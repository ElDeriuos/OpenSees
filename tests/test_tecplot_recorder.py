import os

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


def read_tecplot(path):
    """Minimal reader for the FEBLOCK files the tecplot recorder writes."""
    with open(path) as f:
        text = f.read()
    header, *zones = text.split("ZONE ")
    names = [line.strip().strip('"') for line in header.splitlines()[2:]]
    out = []
    for z in zones:
        head, body = z.split("\n", 1)
        keys = dict(kv.split("=", 1) for kv in head.replace(" ", "").split(",") if "=" in kv)
        n = int(keys["NODES"] if "NODES" in keys else keys["N"])
        e = int(keys["ELEMENTS"] if "ELEMENTS" in keys else keys["E"])
        size = 8 if "BRICK" in head else 4
        tokens = body.split()
        assert len(tokens) == len(names) * n + e * size, "zone size does not match its header"
        values = {name: [float(x) for x in tokens[i * n:(i + 1) * n]] for i, name in enumerate(names)}
        conn = [[int(x) for x in tokens[len(names) * n + k * size:len(names) * n + (k + 1) * size]]
                for k in range(e)]
        assert all(1 <= c <= n for cell in conn for c in cell), "connectivity out of range"
        out.append((head, values, conn))
    return names, out


def _frame_model():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, x, y in [(1, 0, 0), (2, 1, 0), (3, 1, 1), (4, 0, 1), (5, 2, 0), (6, 2, 0), (7, 3, 3)]:
        ops.node(tag, x, y)
    ops.fix(1, 1, 1)
    ops.fix(4, 1, 1)
    ops.fix(6, 1, 1)
    ops.fix(7, 1, 1)  # free-standing node, drawn only as a point
    ops.uniaxialMaterial("Elastic", 1, 1000.0)
    ops.nDMaterial("ElasticIsotropic", 2, 1000.0, 0.25)
    ops.element("Truss", 1, 2, 5, 1.0, 1)
    ops.element("quad", 2, 1, 2, 3, 4, 1.0, "PlaneStress", 2)
    ops.element("zeroLength", 3, 5, 6, "-mat", 1, 1, "-dir", 1, 2)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, 1.0, 0.5)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 0.5)
    ops.algorithm("Linear")
    ops.analysis("Static")


def test_tecplot_recorder_values_match_model(tmp_path):
    path = os.path.join(str(tmp_path), "model.plt")
    _frame_model()
    ops.recorder("tecplot", path, "disp", "eleResponse", "axialForce", "-precision", 12)
    ops.analyze(2)
    disp5 = ops.nodeDisp(5)
    force = ops.eleResponse(1, "axialForce")
    ops.wipe()  # closes the file

    names, zones = read_tecplot(path)
    assert names[:6] == ["X", "Y", "Z", "NodeTag", "EleTag", "EleClassTag"]
    assert "Displacement_1" in names and "Truss_axialForce_1" in names
    assert len(zones) == 2
    head, values, conn = zones[-1]
    assert "SOLUTIONTIME=1.0" in head.replace(" ", "") and "FEQUADRILATERAL" in head

    # 7 nodes as points + truss (2) + quad (4) + zeroLength (2) node copies
    assert len(values["X"]) == 7 + 2 + 4 + 2 and len(conn) == 7 + 3
    for i, tag in enumerate(values["NodeTag"]):
        if tag == 5:
            assert values["Displacement_1"][i] == pytest.approx(disp5[0])
            assert values["Displacement_2"][i] == pytest.approx(disp5[1])
        is_truss = values["EleTag"][i] == 1
        assert values["Truss_axialForce_1"][i] == pytest.approx(force[0] if is_truss else 0.0)

    # degenerate cells: a point repeats one node, the truss is a line
    assert len(set(conn[0])) == 1
    truss = [c for c in conn if values["EleTag"][c[0] - 1] == 1][0]
    assert truss[2:] == [truss[1], truss[1]]


def test_tecplot_recorder_legacy_and_bricks(tmp_path):
    path = os.path.join(str(tmp_path), "solid.plt")
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    coords = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1), (2, 0, 0)]
    for tag, xyz in enumerate(coords, 1):
        ops.node(tag, *xyz)
    for tag in (1, 2, 3, 4):
        ops.fix(tag, 1, 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.element("FourNodeTetrahedron", 2, 2, 9, 3, 6, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(7, 0.0, 0.0, -1.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    ops.recorder("tecplot", path, "disp", "-legacy")
    ops.analyze(1)
    ops.wipe()

    _, zones = read_tecplot(path)
    head, values, conn = zones[0]
    assert "F=FEBLOCK" in head and "ET=BRICK" in head and "STRANDID" not in head
    tet = [c for c in conn if values["EleTag"][c[0] - 1] == 2][0]
    a, b, c, d = tet[0], tet[1], tet[2], tet[4]
    assert tet == [a, b, c, c, d, d, d, d]
