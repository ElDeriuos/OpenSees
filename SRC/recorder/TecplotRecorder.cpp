/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the  **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** ****************************************************************** */

// Description: see TecplotRecorder.h.
//
// recorder tecplot fileName <disp vel accel incrDisp reaction pressure unbalancedLoad mass>
//          <eigen n> <eleResponse "args" ...> <-precision n> <-dT dt> <-rTolDt tol>
//          <-region tag> <-legacy>

#include "TecplotRecorder.h"
#include "PVDRecorder.h"
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Node.h>
#include <NodeIter.h>
#include <Pressure_Constraint.h>
#include <Pressure_ConstraintIter.h>
#include <MeshRegion.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <classTags.h>
#include "PFEMElement/BackgroundDef.h"
#include "PFEMElement/Particle.h"
#include "PFEMElement/ParticleGroup.h"
#include <cstring>
#include <map>
#include <sstream>

// variable kinds besides TecplotRecorder::NodeQuantity
enum {KX = 100, KY, KZ, KNodeTag, KEleTag, KEleClass, KEleResp};

void *OPS_TecplotRecorder()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING: recorder tecplot fileName <disp vel accel incrDisp reaction pressure "
                  "unbalancedLoad mass> <eigen n> <eleResponse \"args\" ...> <-precision n> "
                  "<-dT dt> <-rTolDt tol> <-region tag> <-legacy>\n";
        return 0;
    }
    const char *name = OPS_GetString();

    static const char *quantities[] = {"disp", "vel", "accel", "incrDisp", "reaction",
                                       "unbalancedLoad", "mass", "pressure"};
    std::vector<int> nodeQuantities;
    std::vector<std::vector<std::string> > eleResponses;
    int numEigen = 0, precision = 10, regionTag = -1;
    double dT = 0.0, rTolDt = 0.00001;
    bool legacy = false;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();
        int q = 0;
        while (q < 8 && strcmp(opt, quantities[q]) != 0) q++;
        int one = 1, ok = 0;
        if (q < 8) {
            nodeQuantities.push_back(q);
        } else if (strcmp(opt, "-legacy") == 0) {
            legacy = true;
        } else if (strcmp(opt, "eleResponse") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: recorder tecplot - eleResponse needs its arguments\n";
                return 0;
            }
            // one argument split on spaces, so "section 1 force" works
            std::istringstream ss(OPS_GetString());
            std::vector<std::string> tokens;
            std::string t;
            while (ss >> t) tokens.push_back(t);
            if (!tokens.empty()) eleResponses.push_back(tokens);
        } else if (strcmp(opt, "eigen") == 0) {
            ok = OPS_GetIntInput(&one, &numEigen);
        } else if (strcmp(opt, "-precision") == 0) {
            ok = OPS_GetIntInput(&one, &precision);
        } else if (strcmp(opt, "-region") == 0) {
            ok = OPS_GetIntInput(&one, &regionTag);
        } else if (strcmp(opt, "-dT") == 0) {
            ok = OPS_GetDoubleInput(&one, &dT);
        } else if (strcmp(opt, "-rTolDt") == 0) {
            ok = OPS_GetDoubleInput(&one, &rTolDt);
        } else {
            opserr << "WARNING: recorder tecplot - unknown option " << opt << " ignored\n";
        }
        if (ok < 0) {
            opserr << "WARNING: recorder tecplot - invalid value for " << opt << "\n";
            return 0;
        }
    }
    if (dT < 0) dT = 0;
    if (rTolDt < 0) rTolDt = 0;

    return new TecplotRecorder(name, nodeQuantities, numEigen, eleResponses, precision,
                               dT, rTolDt, regionTag, legacy);
}

TecplotRecorder::TecplotRecorder(const char *name, const std::vector<int> &nq, int neigen,
                                 const std::vector<std::vector<std::string> > &er, int pre,
                                 double dt, double rTolDt, int region, bool leg)
    :Recorder(RECORDER_TAGS_TecplotRecorder), filename(name), nodeQuantities(nq),
     numEigen(neigen), eleResponses(er), precision(pre), regionTag(region), dT(dt),
     relDeltaTTol(rTolDt), nextTime(0.0), legacy(leg), theDomain(0), numDOF(1)
{
}

TecplotRecorder::~TecplotRecorder()
{
    theFile.close();
}

int
TecplotRecorder::writeHeader(const std::vector<Element *> &eles)
{
    static const char *base[] = {"X", "Y", "Z", "NodeTag", "EleTag", "EleClassTag"};
    static const char *names[] = {"Displacement", "Velocity", "Acceleration", "IncrDisplacement",
                                  "Reaction", "UnbalancedLoad", "NodeMass", "Pressure"};
    vars.clear();
    for (int k = 0; k < 6; k++) {
        Var v = {base[k], KX + k, 0, -1, 0};
        vars.push_back(v);
    }

    numDOF = 1;
    NodeIter &theNodes = theDomain->getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0)
        if (theNode->getNumberDOF() > numDOF) numDOF = theNode->getNumberDOF();
    TaggedObjectIter &meshes = OPS_getAllMesh();  // particle-only models have no nodes
    TaggedObject *mesh;
    while ((mesh = meshes()) != 0) {
        ParticleGroup *group = dynamic_cast<ParticleGroup *>(mesh);
        if (group != 0 && group->numParticles() > 0 && (int)group->getParticle(0)->getVel().size() > numDOF)
            numDOF = (int)group->getParticle(0)->getVel().size();
    }

    for (size_t i = 0; i < nodeQuantities.size(); i++) {
        int q = nodeQuantities[i];
        for (int j = 0; j < (q == Pressure ? 1 : numDOF); j++) {
            Var v = {q == Pressure ? names[q] : names[q] + ("_" + std::to_string(j + 1)), q, j, -1, 0};
            vars.push_back(v);
        }
    }
    for (int k = 0; k < numEigen; k++)
        for (int j = 0; j < numDOF; j++) {
            Var v = {"EigenVector" + std::to_string(k + 1) + "_" + std::to_string(j + 1), Eigen, j, k, 0};
            vars.push_back(v);
        }

    // element responses are sized from the first element of each class, as in PVDRecorder.
    // shortcut: a class that first appears after this record gets no response columns.
    std::map<int, Element *> firstOfClass;
    for (size_t i = 0; i < eles.size(); i++) firstOfClass.insert(std::make_pair(eles[i]->getClassTag(), eles[i]));
    for (size_t r = 0; r < eleResponses.size(); r++) {
        std::vector<const char *> argv;
        std::string suffix;
        for (size_t t = 0; t < eleResponses[r].size(); t++) {
            argv.push_back(eleResponses[r][t].c_str());
            suffix += "_" + eleResponses[r][t];
        }
        for (std::map<int, Element *>::iterator it = firstOfClass.begin(); it != firstOfClass.end(); it++) {
            const Vector *data = theDomain->getElementResponse(it->second->getTag(), &argv[0], (int)argv.size());
            if (data == 0) continue;
            for (int k = 0; k < data->Size(); k++) {
                Var v = {it->second->getClassType() + suffix + "_" + std::to_string(k + 1), KEleResp, k, (int)r, it->first};
                vars.push_back(v);
            }
        }
    }

    theFile.open(filename.c_str(), std::ios::trunc | std::ios::out);
    if (theFile.fail()) {
        opserr << "WARNING: recorder tecplot - failed to open " << filename.c_str() << "\n";
        return -1;
    }
    theFile.precision(precision);
    theFile << std::scientific << "TITLE = \"OpenSees\"\nVARIABLES =";
    // one name per line keeps every line under Tecplot 8's header line limit
    for (size_t i = 0; i < vars.size(); i++) theFile << "\n\"" << vars[i].name << "\"";
    theFile << "\n";
    return 0;
}

double
TecplotRecorder::value(const Slot &s, const Var &v, const std::vector<std::vector<double> > &resp)
{
    switch (v.kind) {
    case KX: case KY: case KZ: {
        int j = v.kind - KX;
        if (s.particle != 0) {
            const VDouble &c = s.particle->getCrds();
            return j < (int)c.size() ? c[j] : 0.0;
        }
        const Vector &c = s.node->getCrds();
        return j < c.Size() ? c(j) : 0.0;
    }
    case KNodeTag: return s.particle != 0 ? (double)s.particle->getTag() : s.node->getTag();
    case KEleTag: return s.ele != 0 ? s.ele->getTag() : 0;
    case KEleClass: return s.ele != 0 ? s.ele->getClassTag() : 0;
    case KEleResp: {
        if (s.ele == 0 || s.ele->getClassTag() != v.classTag) return 0.0;
        const std::vector<double> &r = resp[s.eleIndex * eleResponses.size() + v.resp];
        return v.comp < (int)r.size() ? r[v.comp] : 0.0;
    }
    }

    if (s.particle != 0) {  // particles carry only velocity and pressure
        if (v.kind == Vel) {
            const VDouble &u = s.particle->getVel();
            return v.comp < (int)u.size() ? u[v.comp] : 0.0;
        }
        return v.kind == Pressure ? s.particle->getPressure() : 0.0;
    }

    Node *nd = s.node;
    const Vector *u = 0;
    switch (v.kind) {
    case Disp: u = &nd->getTrialDisp(); break;
    case Vel: u = &nd->getTrialVel(); break;
    case Accel: u = &nd->getTrialAccel(); break;
    case IncrDisp: u = &nd->getIncrDisp(); break;
    case Reaction: u = &nd->getReaction(); break;
    case Unbalanced: u = &nd->getUnbalancedLoad(); break;
    case Mass: {
        const Matrix &m = nd->getMass();
        return v.comp < m.noRows() ? m(v.comp, v.comp) : 0.0;
    }
    case Pressure: {
        Pressure_Constraint *pc = theDomain->getPressure_Constraint(nd->getTag());
        return pc != 0 ? pc->getPressure() : 0.0;
    }
    case Eigen: {
        // Node::getEigenvectors() exits the program when no eigen analysis has run
        if (v.resp >= theDomain->getEigenvalues().Size()) return 0.0;
        const Matrix &m = nd->getEigenvectors();
        return (v.comp < m.noRows() && v.resp < m.noCols()) ? m(v.comp, v.resp) : 0.0;
    }
    }
    return (u != 0 && v.comp < u->Size()) ? (*u)(v.comp) : 0.0;
}

int
TecplotRecorder::record(int commitTag, double timeStamp)
{
    if (dT > 0.0) {
        if (timeStamp - nextTime < -dT * relDeltaTTol) return 0;
        nextTime = timeStamp + dT;
    }
    if (theDomain == 0) return -1;

    // elements to record: a region, else the whole domain
    std::vector<Element *> eles;
    MeshRegion *region = regionTag >= 0 ? theDomain->getRegion(regionTag) : 0;
    if (regionTag >= 0 && region == 0 && warned.insert(-1).second)
        opserr << "WARNING: recorder tecplot - region " << regionTag << " not found, recording the whole domain\n";
    Element *theEle;
    if (region != 0) {
        const ID &tags = region->getElements();
        for (int i = 0; i < tags.Size(); i++)
            if ((theEle = theDomain->getElement(tags(i))) != 0) eles.push_back(theEle);
    } else {
        ElementIter &theEles = theDomain->getElements();
        while ((theEle = theEles()) != 0) eles.push_back(theEle);
    }

    if (!theFile.is_open() && this->writeHeader(eles) < 0) return -1;
    // shortcut: like pvd, asking for eigen before any eigen analysis exits inside
    // Domain::getEigenvalues(); a public Domain "has eigenvalues" accessor would avoid it.
    if (numEigen > 0 && numEigen > theDomain->getEigenvalues().Size() && warned.insert(-2).second)
        opserr << "WARNING: recorder tecplot - fewer than " << numEigen << " eigenvectors available, writing zeros\n";

    // every node, particle and element becomes a cell over its own slots (node copies)
    struct Cell {int first, count, type;};
    std::vector<Slot> slots;
    std::vector<Cell> cells;

    // nodes as point cells, like part 0 of PVDRecorder; PFEM pressure nodes are internal
    std::set<int> skip, regionNodes;
    Pressure_ConstraintIter &thePCs = theDomain->getPCs();
    Pressure_Constraint *thePC;
    while ((thePC = thePCs()) != 0)
        if (thePC->getPressureNode() != 0) skip.insert(thePC->getPressureNode()->getTag());
    for (size_t i = 0; region != 0 && i < eles.size(); i++) {
        const ID &n = eles[i]->getExternalNodes();
        for (int j = 0; j < n.Size(); j++) regionNodes.insert(n(j));
    }
    NodeIter &theNodes = theDomain->getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0) {
        int tag = theNode->getTag();
        if (skip.count(tag) || (region != 0 && !regionNodes.count(tag))) continue;
        Slot s = {theNode, 0, 0, -1};
        Cell c = {(int)slots.size(), 1, PVDRecorder::VTK_VERTEX};
        slots.push_back(s);
        cells.push_back(c);
    }

    // PFEM particles (not part of regions)
    TaggedObjectIter &meshes = OPS_getAllMesh();
    TaggedObject *mesh;
    while (region == 0 && (mesh = meshes()) != 0) {
        ParticleGroup *group = dynamic_cast<ParticleGroup *>(mesh);
        for (int j = 0; group != 0 && j < group->numParticles(); j++) {
            Slot s = {0, group->getParticle(j), 0, -1};
            Cell c = {(int)slots.size(), 1, PVDRecorder::VTK_VERTEX};
            slots.push_back(s);
            cells.push_back(c);
        }
    }

    // elements by their corner nodes; responses copied now (Domain reuses its buffer)
    std::set<int> respClasses;
    for (size_t i = 0; i < vars.size(); i++)
        if (vars[i].kind == KEleResp) respClasses.insert(vars[i].classTag);
    std::vector<std::vector<double> > resp(eles.size() * eleResponses.size());
    bool brick = false;
    for (size_t i = 0; i < eles.size(); i++) {
        Element *e = eles[i];
        int ctag = e->getClassTag();
        int type = PVDRecorder::getVTKType(ctag);
        if (type == 0) {
            if (warned.insert(ctag).second)
                opserr << "WARNING: recorder tecplot - no shape for element class " << e->getClassType() << ", not drawn\n";
            continue;
        }
        int num = 0;
        switch (type) {
        case PVDRecorder::VTK_LINE: num = 2; break;
        case PVDRecorder::VTK_TRIANGLE: case PVDRecorder::VTK_QUADRATIC_TRIANGLE: num = 3; break;
        case PVDRecorder::VTK_QUAD: case PVDRecorder::VTK_QUADRATIC_QUAD:
        case PVDRecorder::VTK_TETRA: case PVDRecorder::VTK_QUADRATIC_TETRA: num = 4; break;
        case PVDRecorder::VTK_HEXAHEDRON: case PVDRecorder::VTK_QUADRATIC_HEXAHEDRON: num = 8; break;
        default: num = 4;  // vertex types (zeroLength, joints, 9-node quads ...): first nodes as point/line/tri/quad
        }
        // PFEM elements interleave pressure nodes (same list as PVDRecorder::savePart)
        int stride = (ctag == ELE_TAG_PFEMElement2D || ctag == ELE_TAG_PFEMElement2DCompressible ||
                      ctag == ELE_TAG_PFEMElement2DBubble || ctag == ELE_TAG_PFEMElement2Dmini ||
                      ctag == ELE_TAG_MINI || ctag == ELE_TAG_PFEMElement2DQuasi ||
                      ctag == ELE_TAG_PFEMElement3DBubble) ? 2 : 1;
        const ID &n = e->getExternalNodes();
        std::vector<Node *> corners;
        for (int j = 0; j < num && j * stride < n.Size(); j++) {
            Node *nd = theDomain->getNode(n(j * stride));
            if (nd == 0) break;
            corners.push_back(nd);
        }
        bool vertexType = type == PVDRecorder::VTK_POLY_VERTEX || type == PVDRecorder::VTK_VERTEX;
        if (corners.empty() || ((int)corners.size() < num && !vertexType)) continue;
        if (type == PVDRecorder::VTK_TETRA || type == PVDRecorder::VTK_QUADRATIC_TETRA ||
            type == PVDRecorder::VTK_HEXAHEDRON || type == PVDRecorder::VTK_QUADRATIC_HEXAHEDRON) brick = true;

        Cell c = {(int)slots.size(), (int)corners.size(), type};
        for (size_t j = 0; j < corners.size(); j++) {
            Slot s = {corners[j], 0, e, (int)i};
            slots.push_back(s);
        }
        cells.push_back(c);

        for (size_t r = 0; respClasses.count(ctag) && r < eleResponses.size(); r++) {
            std::vector<const char *> argv;
            for (size_t t = 0; t < eleResponses[r].size(); t++) argv.push_back(eleResponses[r][t].c_str());
            const Vector *data = theDomain->getElementResponse(e->getTag(), &argv[0], (int)argv.size());
            for (int k = 0; data != 0 && k < data->Size(); k++) resp[i * eleResponses.size() + r].push_back((*data)(k));
        }
    }
    if (cells.empty()) return 0;

    // zone header
    theFile << "ZONE T=\"t = " << timeStamp << "\", ";
    if (legacy) {
        theFile << "N=" << slots.size() << ", E=" << cells.size() << ", F=FEBLOCK, ET="
                << (brick ? "BRICK" : "QUADRILATERAL") << "\n";
    } else {
        theFile << "STRANDID=1, SOLUTIONTIME=" << timeStamp << ", ZONETYPE="
                << (brick ? "FEBRICK" : "FEQUADRILATERAL") << ", NODES=" << slots.size()
                << ", ELEMENTS=" << cells.size() << ", DATAPACKING=BLOCK, DT=(";
        for (size_t i = 0; i < vars.size(); i++) theFile << "DOUBLE ";
        theFile << ")\n";
    }

    // one block per variable, 10 values per line
    for (size_t k = 0; k < vars.size(); k++)
        for (size_t i = 0; i < slots.size(); i++)
            theFile << this->value(slots[i], vars[k], resp) << ((i + 1) % 10 == 0 || i + 1 == slots.size() ? '\n' : ' ');

    // connectivity (1-based): short cells repeat their last node, Tecplot's degenerate
    // elements; a tetrahedron in a brick zone uses Tecplot's tet pattern instead
    static const int tet[8] = {0, 1, 2, 2, 3, 3, 3, 3};
    int size = brick ? 8 : 4;
    for (size_t i = 0; i < cells.size(); i++) {
        bool isTet = cells[i].type == PVDRecorder::VTK_TETRA || cells[i].type == PVDRecorder::VTK_QUADRATIC_TETRA;
        for (int j = 0; j < size; j++) {
            int local = (brick && isTet) ? tet[j] : (j < cells[i].count ? j : cells[i].count - 1);
            theFile << cells[i].first + local + 1 << (j + 1 < size ? ' ' : '\n');
        }
    }
    theFile.flush();
    return 0;
}

int
TecplotRecorder::restart(void)
{
    theFile.close();  // the next record starts the file again
    return 0;
}

int
TecplotRecorder::domainChanged(void)
{
    return 0;  // zones are rebuilt from the domain at every record
}

int
TecplotRecorder::setDomain(Domain &domain)
{
    theDomain = &domain;
    return 0;
}

int
TecplotRecorder::flush(void)
{
    if (theFile.is_open()) theFile.flush();
    return 0;
}

int
TecplotRecorder::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
TecplotRecorder::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}
