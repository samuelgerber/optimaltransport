#ifndef TRANSPORTLPSOLVER_H
#define TRANSPORTLPSOLVER_H


#include "MultiscaleTransport.h"
#include "LPSolver.h"


template <typename TPrecision>
class TransportLPSolver {


  public:

    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;
    typedef typename TransportPlan<TPrecision>::Path Path;
    typedef LPSolver::Status Status;


    static void createLP(TransportPlan<TPrecision> *sol, LPSolver *solver){
      int ns = sol->source->getNodes().size();
      int nt = sol->target->getNodes().size();

      solver->createLP(ns, nt);

      //rows - constraints:
      // M of type sum of weights going from node A tree1 to nodes in tree2 =
      // weight of node A in tree1
      // N of type sum of weights from nodes in tree1 going to node A in tree2
      // = weight of node A in tree2
      int nConstraints = ns + nt;

      solver->addRows( nConstraints );

      //columns - coefficents for each combination of nodes from tree1 to
      //tree2 = M*N
      solver->addCols( sol->getNumberOfPaths() );

      //set constraints bounds
      //The weights of each node in tree1
      for(TransportNodeVectorCIterator it = sol->source->getNodes().begin(); it !=
          sol->source->getNodes().end(); ++it){
        TransportNode<TPrecision> *n = *it;
        double w = n->getMass();
        solver->setRowBounds(n->getID(), w);
      }

      int offset= sol->source->getNodes().size();
      //the weight of each node in tree2
      for(TransportNodeVectorCIterator it = sol->target->getNodes().begin(); it !=
          sol->target->getNodes().end(); ++it){
        TransportNode<TPrecision> *n = *it;
        double w = n->getMass();
        solver->setRowBounds( offset + n->getID(), -w );
      }




      //Setup constraint matrix based on the following set of constraints:
      //1. the sum of the N weights leaving from A has to equal the weight in the
      //receiver in B
      //2. the sum of the M weights ending up in node i of tree2 has to equal the
      //weight of node i in tree2
      //The constraint bounds above are the receiving and sending amount of
      //mass. Each path needs to add a 1 in the corresponding rows (from and
      //to row)

      for(sol->pathIteratorBegin(); !sol->pathIteratorIsAtEnd();
          sol->pathIteratorNext() ){
        Path &path = sol->pathIteratorCurrent();
        TransportNode<TPrecision> *from = path.from;
        TransportNode<TPrecision> *to = path.to;
        solver->setColBounds(path.index);
        solver->setCoefficent(path.index, path.cost );
        solver->setColConstraints(path.index, from->getID(), offset + to->getID() );

      }

   }




   static void storeLP(TransportPlan<TPrecision> *sol, LPSolver *solver, TPrecision p){
     sol->cost = pow( solver->getObjectiveValue(), 1.0/p );


     std::cout << "cost: " << sol->cost << std::endl;
     //Store solution
     double sumw = 0;
     int nNonZero =0;
     for(sol->pathIteratorBegin(); !sol->pathIteratorIsAtEnd();
         sol->pathIteratorNext() ){
       Path &path = sol->pathIteratorCurrent();
       path.w = solver->getColPrimal(path.index);
       sumw += path.w;
       nNonZero += (path.w > 0);
     }

     std::cout << "nonzeros: "      << nNonZero << std::endl;
     std::cout << "solution sumw: " << sumw     << std::endl << std::endl;


     //Potential of from nodes
     for(TransportNodeVectorCIterator it = sol->source->getNodes().begin(); it !=
         sol->source->getNodes().end(); ++it ){
       TransportNode<TPrecision> *node = *it;
       TPrecision pi = solver->getRowDual( node->getID() );
       node->setPotential(pi);
     }


     //Potential of to nodes
     int offset = sol->source->getNodes().size();
     for(TransportNodeVectorCIterator it = sol->target->getNodes().begin(); it !=
         sol->target->getNodes().end(); ++it){
       TransportNode<TPrecision> *node = *it;
       TPrecision pi = solver->getRowDual( offset + node->getID() );
       node->setPotential(pi);
     }

   };







   static void setPotentials(TransportNodeVector &sourceNodes, TransportNodeVector
        &targetNodes, LPSolver *solver){

      for(TransportNodeVectorIterator sIt = sourceNodes.begin(); sIt !=
          sourceNodes.end(); ++sIt){
        TransportNode<TPrecision> *node = *sIt;
        TPrecision pi = solver->getRowDual(node->getID());
        node->setPotential(pi);
        node->resetPi(pi, pi);
      }

      int offset = sourceNodes.size();
      for(TransportNodeVectorIterator tIt = targetNodes.begin(); tIt !=
          targetNodes.end(); ++tIt){
        TransportNode<TPrecision> *node = *tIt;
        TPrecision pi = solver->getRowDual( offset + node->getID() );
        node->setPotential(pi);
        node->resetPi(pi, pi);
      }

    };




    static void setupBasis(LPSolver *solver, std::vector<Status> &colStatus,
        std::vector<Status> &rowStatus){

      //int nBasic = 0;
      for(int i=0; i<colStatus.size(); i++){
        solver->setColStatus(i, colStatus[i]);
        //nBasic += colStatus[i] == LPSolver::BASIC;
      }

      for(int i=0; i<rowStatus.size(); i++){
        solver->setRowStatus(i, rowStatus[i]);
        //nBasic += rowStatus[i] == LPSolver::BASIC;
      }

     // std::cout << "#basic: " << nBasic << std::endl;

    };




};


#endif


