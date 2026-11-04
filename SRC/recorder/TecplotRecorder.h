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

// Description: TecplotRecorder writes the model and its results to one
// ASCII Tecplot file (e.g. model.plt), one FE zone per recorded step.
// Every node, particle and element becomes a (possibly degenerate)
// quadrilateral or brick with its own copies of the nodes. That keeps
// element results piecewise constant without cell-centered data, so the
// same layout loads in current Tecplot 360 (time strands: animation and
// probe-over-time plots) and, with -legacy, in Tecplot 8/9 (Animate Zones).

#ifndef TecplotRecorder_h
#define TecplotRecorder_h

#include <Recorder.h>
#include <fstream>
#include <set>
#include <string>
#include <vector>

class Domain;
class Node;
class Element;
class Particle;

class TecplotRecorder : public Recorder
{
public:
    enum NodeQuantity {Disp, Vel, Accel, IncrDisp, Reaction, Unbalanced, Mass, Pressure, Eigen};

    TecplotRecorder(const char *filename, const std::vector<int> &nodeQuantities, int numEigen,
                    const std::vector<std::vector<std::string> > &eleResponses, int precision,
                    double dT, double rTolDt, int regionTag, bool legacy);
    ~TecplotRecorder();

    int record(int commitTag, double timeStamp);
    int restart(void);
    int domainChanged(void);
    int setDomain(Domain &theDomain);
    int flush(void);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

private:
    struct Var {std::string name; int kind; int comp; int resp; int classTag;};
    struct Slot {Node *node; Particle *particle; Element *ele; int eleIndex;};

    int writeHeader(const std::vector<Element *> &eles);
    double value(const Slot &s, const Var &v, const std::vector<std::vector<double> > &resp);

    std::string filename;
    std::vector<int> nodeQuantities;
    int numEigen;
    std::vector<std::vector<std::string> > eleResponses;
    int precision, regionTag;
    double dT, relDeltaTTol, nextTime;
    bool legacy;

    Domain *theDomain;
    std::ofstream theFile;
    std::vector<Var> vars;      // fixed at the first record: Tecplot needs one variable list per file
    int numDOF;                 // components written per nodal quantity
    std::set<int> warned;       // classes / issues already reported
};

#endif
