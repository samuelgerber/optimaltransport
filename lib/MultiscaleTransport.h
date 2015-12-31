#ifndef MULTISCALETRANSPORT_H
#define MULTISCALETRANSPORT_H

#include <vector>
#include <list>
#include <map>
#include <ctime>
//#include <boost/unordered_map.hpp>
//#include <boost/unordered_set.hpp>

#include "Eigen/Dense"


template <typename TPrecision>
class TransportNode{
  
  public:
    typedef typename std::vector< TransportNode<TPrecision>*  >  TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
  
  private:
  
    TPrecision piMax;
    TPrecision piMin;
    TPrecision potential;

    TPrecision mass;

    TransportNodeVector kids;
    TransportNode<TPrecision> *parent; 

    TPrecision radius;

    int id;

    int scale;


  public:

    TransportNode(int nodeID, int sca):id(nodeID), scale(sca){
      piMin = std::numeric_limits<TPrecision>::max();
      piMax = -piMin;
      mass = -1;
      parent = NULL;
      radius = -1;
    };

    virtual ~TransportNode(){
    };


    TPrecision getMass() const{
      return mass;
    };


    void setMass(TPrecision m){
      //std::cout << m << std::endl;
      mass = m;
    };


    void addChild(TransportNode *node){
      kids.push_back(node);
    };


    TransportNodeVector getChildren() const{
      return kids;
    };


    void setParent(TransportNode<TPrecision> *p){
      parent = p;
    };


    TransportNode<TPrecision> *getParent() const{
      return parent;
    };


    virtual TPrecision getNodeRadius() const = 0;
    virtual TPrecision getLocalNodeRadius() const = 0;
    virtual std::vector<int> &getPoints() = 0;

    virtual TPrecision getTransportCostRadius(double p){

      //NOTE: will *not* update if children are added after calling getRadius
      if(radius < 0){
        radius = 0;
        for(TransportNodeVectorIterator it = kids.begin(); it != kids.end(); ++it){
          TPrecision d = getTransportCost(*it, p);
          if(d > radius){
            radius = d;
          }
        }
      }

      return radius;
    };

    
    TPrecision getAverageTransportCost(const TransportNode<TPrecision> *other, double p){
      TPrecision sum = 0;
      for(int i=0; i<kids.size(); i++){
        for(int j=0; j<other->kids.size(); j++){
          sum += kids[i]->getTransportCost(other->kids[j], p);
        }
      }

      return sum / ( kids.size() * other->kids.size() );
    
    };

    virtual TPrecision getTransportCost(const TransportNode<TPrecision> *other, double p) const= 0;
    //virtual TPrecision getDeltaTransportCost(TransportNode<TPrecision> *other, double p) const= 0;
    
    /*
    virtual bool operator == (const TransportNode<TPrecision> &other) const = 0;
    virtual bool operator <  (const TransportNode<TPrecision> &other) const = 0;
    virtual bool operator >  (const TransportNode<TPrecision> &other) const = 0;
    */
    //virtual size_t hashKey() const = 0; 



    void resetPi(double piMi = std::numeric_limits<TPrecision>::max(), double
        piMa = - std::numeric_limits<TPrecision>::max() ){
      piMin = piMi;
      piMax = piMa;
    };

    void setPiMin(TPrecision pi){
      if(piMin > pi){
        piMin = pi;
      }
    };

    void setPiMax(TPrecision pi){
      if(piMax < pi){
        piMax = pi;
      }
    };

    TPrecision getPiMin(){
      return piMin;
    };

    TPrecision getPiMax(){
      return piMax;
    };

    void setPotential(TPrecision p){
      potential = p;
    };

    TPrecision getPotential() const{
      return potential;
    };

    int getID() const{
      return id;
    };

    int getScale() const{
      return scale;
    };


};





template <typename TPrecision>
class TransportNodeEqual{
  public:
    bool operator()(const TransportNode<TPrecision> *s1, const
        TransportNode<TPrecision> *s2) const{
      return s1->operator == (*s2);
    };
};




template <typename TPrecision>
class TransportNodeLess{
  public:
    bool operator()(const TransportNode<TPrecision> *s1, const
        TransportNode<TPrecision> *s2) const{
      return s1->operator < (*s2);
    };
};




template <typename TPrecision>
class TransportNodeHash{
  public:
    bool operator()(const TransportNode<TPrecision> *s1) const{
      return s1->hashKey();
    };
};






template <typename TPrecision>
class MultiscaleTransportLevel{
 
  public:
  
    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;
  

  private:
  
    int scale;
    TransportNodeVector nodes;
 
    MultiscaleTransportLevel<TPrecision> *parent;

  public:
  
    MultiscaleTransportLevel(int s, MultiscaleTransportLevel<TPrecision> *parentLevel) : scale(s), parent(parentLevel){
    };


    virtual ~MultiscaleTransportLevel(){
      for(TransportNodeVectorIterator it = nodes.begin(); it != nodes.end(); ++it){
        delete *it;
      }
      nodes.clear();
    };


    void addNode(TransportNode<TPrecision> *node){
      nodes.push_back(node);
    };


    TransportNodeVector &getNodes(){
      return nodes;
    };


    virtual TransportNodeVector getNeighborhood(TransportNode<TPrecision> *node,
        TPrecision eps) const = 0;

    TPrecision getMaximalRadius(){
      TPrecision radius = 0;
      TransportNodeVector &nodes = getNodes();
      for(TransportNodeVectorCIterator it = nodes.begin(); it != nodes.end();
          ++it){
        if(radius <  (*it)->getRadius() ){
          radius = (*it)->getRadius();
        }
      }
      return radius;

    };



    int getScale(){
      return scale;
    };

   

    MultiscaleTransportLevel<TPrecision> *getRootLevel(){
      MultiscaleTransportLevel<TPrecision> *current = this;
      MultiscaleTransportLevel<TPrecision> *parent = current->getParent();
      while(parent!=NULL){
        current=parent;
        parent = current->getParent();
      }

      return current;
    };



    MultiscaleTransportLevel<TPrecision> *getParent(){
      return parent;
    };

};














template <typename TPrecision>
class TransportPlan{
 


  public:


    struct Path{
      TransportNode<TPrecision> *from;
      TransportNode<TPrecision> *to;

      Path() {
        from = NULL;
        to = NULL;
        index= -1;
        cost = -1;
        w = 0;
      };


      Path(TransportNode<TPrecision> *f, TransportNode<TPrecision> *t) 
        :from(f),to(t){
          index= -1;
          cost = -1;
          w = 0;
        };


      bool operator < (const Path &s1) const{
        if(from->getID() < s1.from->getID() ){
          return true;
        }
        else if(from->getID() > s1.from->getID() ){
          return false;
        }
        else if(to->getID() < s1.to->getID() ){
          return true;
        }
        return false;
      };


      bool operator > (const Path &s1) const {
        if(from->getID() > s1.from->getID() ){
          return true;
        }
        else if(from->getID() < s1.from->getID() ){
          return false;
        }
        else if(to->getID() > s1.to->getID() ){
          return true;
        }
        return false;
      };

      bool operator == (const Path &s1) const{
        return from->getID() == s1.from->getID() && to->getID() == s1.to->getID();
      };


      int index;
      TPrecision cost;
      TPrecision w;

    };




  private:

    std::vector< std::map< int, Path> > paths;
    typename std::vector< std::map<int, Path> >::iterator pathIteratorOuter;
    typename std::map< int, Path>::iterator pathIteratorInner;

    std::vector<int> toPathCounts;
    int pathCounter;
  


  public:


    TransportPlan(MultiscaleTransportLevel<TPrecision> *s,
        MultiscaleTransportLevel<TPrecision> *t) : source(s), target(t){
      
      paths.resize( source->getNodes().size() );
      toPathCounts.resize( target->getNodes().size(), 0 );
      timeSolve = 0;
      timeRefine = 0;
      timePropagate = 0;
      cost = std::numeric_limits<TPrecision>::max();
      optimizationStatus = -1;
      //nTotalPaths = -1;
      pathCounter = 0;
    
    };    
   
/* 
    TransportPlan(int nFrom, int nTo) : source(NULL), target(NULL){
      
      paths.resize( nFrom );
      toPathCounts.resize( nTo, 0);
      timeLP = 0;
      timeNeighborhood = 0;
      timePropagate = 0;
      cost = std::numeric_limits<TPrecision>::max();
      optimizationStatus = -1;
      nTotalPaths = -1;
      pathCounter = 0;
    };

*/

    virtual ~TransportPlan(){
    };


    MultiscaleTransportLevel<TPrecision> *source;
    MultiscaleTransportLevel<TPrecision> *target;


    double cost;
    int optimizationStatus;
    //int nTotalPaths;
    clock_t timeSolve;
    clock_t timePropagate;
    clock_t timeRefine;


    //For sinkhorn transport 
    Eigen::MatrixXd leftScaling;


    bool hasPath(const Path &p){
       return hasPath( p.from->getID(), p.to->getID() );
    };


    bool hasPath(int from, int to){
      return paths[from].find(to) != paths[from].end();
    };


    Path &getPath(TransportNode<TPrecision> *from, TransportNode<TPrecision> *to){
       return getPath( from->getID(), to->getID() );
    };


    Path &getPath(int from, int to){
      return paths[from][to];
    };
    

    int getNumberOfToPaths(int from){
      return paths[from].size();
    };
    int getNumberOfFromPaths(int to){
      return toPathCounts[to];
    };

    int addPath(Path p){
      if( !hasPath(p) ){
        p.index = pathCounter;
        ++pathCounter;
        paths[p.from->getID()][p.to->getID()] = p;
        toPathCounts[p.to->getID()] += 1;
        return p.index;
      }
      else{
        Path &path = getPath(p.from, p.to);
        path.w = std::max(path.w, p.w);
        return path.index;
      }
    };


    int getPathIndex(const Path &p){
      int from = p.from->getID();
      int to = p.to->getID();
      
      typename std::map<int, Path>::const_iterator it = paths[from].find(to);

      if(it == paths[from].end()){
        return -1;
      }
      return it->second.index;
    };


    void pathIteratorBegin(){
      
      pathIteratorOuter = paths.begin();
      if(pathIteratorOuter != paths.end() ){

        pathIteratorInner = (*pathIteratorOuter).begin();
        
        while( pathIteratorInner == (*pathIteratorOuter).end() ){
          ++pathIteratorOuter;
          if( pathIteratorOuter != paths.end()){
            pathIteratorInner = (*pathIteratorOuter).begin();
          }
          else{
            break;
          }
        }
      }

    }; 


    bool pathIteratorIsAtEnd(){
      return pathIteratorOuter == paths.end();
    };


    Path &pathIteratorCurrent(){
      return pathIteratorInner->second;
    };




    void pathIteratorNext(bool erase=false){
      if(erase){
        (*pathIteratorOuter).erase(pathIteratorInner++);
        pathCounter--;
      }
      else{
        ++pathIteratorInner;
      }
      while( pathIteratorInner == (*pathIteratorOuter).end() ){
        ++pathIteratorOuter;
        if( pathIteratorOuter != paths.end()){
          pathIteratorInner = (*pathIteratorOuter).begin();
        }
        else{
          break;
        }
      }      
    };



    int getNumberOfPaths(){
      return pathCounter;
    };




    //Compute a multicsale path cost
    std::vector<TPrecision> getMultiscaleTransportCost(int p){
      std::vector<TPrecision> costs( std::max(source->getScale(), target->getScale())+1, 0 );

      for( this->pathIteratorBegin(); !this->pathIteratorIsAtEnd();
           this->pathIteratorNext() ){
        Path &path = this->pathIteratorCurrent();
       
        if(path.w > 0 ){
          typename std::vector< TransportNode<TPrecision>* > fNodes;
          TransportNode<TPrecision> *from = path.from;
          while(from != NULL){
            fNodes.push_back(from);
            from = from->getParent();
          }  
          typename std::vector< TransportNode<TPrecision>* > tNodes;
          TransportNode<TPrecision> *to = path.to;
          while(to != NULL){
            tNodes.push_back(to);
            to = to->getParent();
          }

         
          typename std::vector< TransportNode<TPrecision>* >::reverse_iterator fIt = fNodes.rbegin();
          typename std::vector< TransportNode<TPrecision>* >::reverse_iterator tIt = tNodes.rbegin();
          from = *fIt;
          to = *tIt;

          int index = 0;
          while(tIt != tNodes.rend() || fIt != fNodes.rend()){
            costs[index] = costs[index] + from->getTransportCost(to, p) * path.w ;
            ++index;

            if(tIt != tNodes.rend()){
              ++tIt;
            }
            if(tIt != tNodes.rend() ){
              to = *tIt;
            }
               
            if(fIt != fNodes.rend()){
              ++fIt;
            }
            if(fIt != fNodes.rend() ){
              from = *fIt;
            }
       
          }
          
        }

      }

      return costs;
    };
    

    TransportPlan<TPrecision> *createCopy(){
      TransportPlan<TPrecision> *res = new
        TransportPlan<TPrecision>(source, target);

      res->timePropagate = this->timePropagate;
      res->timeSolve = this->timeSolve;
      res->timeRefine = this->timeRefine;
      res->paths = this->paths;
      res->pathCounter = this->pathCounter;
      res->toPathCounts = this->toPathCounts;
      return res;
    };

    void copyFrom(TransportPlan<TPrecision> *res){

      this->timePropagate = res->timePropagate;
      this->timeSolve = res->timeSolve;
      this->timeRefine = res->timeRefine;
      this->paths = res->paths;
      this->pathCounter = res->pathCounter;
      this->toPathCounts = res->toPathCounts;
      return res;
    };

};
















template <typename TPrecision>
class TransportPlanSolutions{
  
  private:  
  
    typedef typename  std::list< TransportPlan<TPrecision> *>::iterator TPIterator;
    typedef typename TransportPlan<TPrecision>::Path Path; 
    
    
    std::list< TransportPlan<TPrecision> *> alternatives;
    TransportPlan<TPrecision> *sol;
 


  public:    
    
    


    TransportPlanSolutions(MultiscaleTransportLevel<TPrecision> *s,
        MultiscaleTransportLevel<TPrecision> *t) {
      sol = new TransportPlan<TPrecision>(s, t);
    }; 


    virtual ~TransportPlanSolutions(){

      for(TPIterator it = alternatives.begin(); it !=alternatives.end(); ++it){
        delete *it;
      }
      alternatives.clear();
    }; 


    void setPrimarySolution(TransportPlan<TPrecision> *primary){
      sol = primary;
    };
    
    TransportPlan<TPrecision> *getPrimarySolution(){
      return sol;
    };


    std::list< TransportPlan<TPrecision> * > &getAlternativeSolutions(){
      return alternatives;
    };


    void addAlternativeSolution( TransportPlan<TPrecision> *a){
      alternatives.push_back(a);
    };



    TransportPlan<TPrecision> *getCombinedPaths(){

      TransportPlan<TPrecision> *res = sol->createCopy(); 

      for(TPIterator it = alternatives.begin(); it !=alternatives.end(); ++it){
        TransportPlan<TPrecision> *a= *it;
        for( a->pathIteratorBegin(); !a->pathIteratorIsAtEnd();
                 a->pathIteratorNext() ){
          Path &path = a->pathIteratorCurrent();
          res->addPath(path);
        }
      }

      return res;
    };

};







template <typename TPrecision>
class MultiscaleTransport{

  public:

    
    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;


    typedef typename TransportPlan<TPrecision>::Path Path; 


    virtual ~MultiscaleTransport(){};


    std::vector< TransportPlan<TPrecision> * > solve(std::vector<
        MultiscaleTransportLevel<TPrecision>* > &aLevels,
        std::vector<MultiscaleTransportLevel<TPrecision>* > &bLevels, 
        double p = 1, int nScales1 = -1, int nScales2 = -1, bool matchStartLevel =
        false){

      //compute distances at each tree level
      typename std::vector< TransportPlan<TPrecision> *> solutions; 

      if(nScales1 < 0 || nScales1 >= (int) aLevels.size()){
        nScales1 = aLevels.size()-1;
      }
      
      if(nScales2 < 0 || nScales2 >= (int) bLevels.size()){
        nScales2 = bLevels.size()-1;
      }

      normalize(aLevels);
      normalize(bLevels);

      typename std::vector< MultiscaleTransportLevel<TPrecision> * >::iterator itA =
        aLevels.begin();
      rootSource = *itA;

      typename std::vector< MultiscaleTransportLevel<TPrecision> * >::iterator itB =
        bLevels.begin();
      rootTarget = *itB;




      int scaleStart1 = aLevels.size()-1-nScales1;
      int scaleStart2 = bLevels.size()-1-nScales2;
      std::advance(itA, scaleStart1);
      std::advance(itB, scaleStart2);

      if(matchStartLevel){

        TPrecision rA = getMeanRadius(*itA, p);
        TPrecision rB = getMeanRadius(*itB, p);
        TPrecision delta = rA-rB;
        if(delta > 0){
          ++itA;
          while(itA != aLevels.end() ){
            rA = getMeanRadius(*itA, p);
            TPrecision tmp = rA-rB;
            if( fabs(tmp) > fabs(delta) ){
              --itA;
              break;
            }
            else{
              delta = tmp;
            }
            ++itA;
          }

          if(itA == aLevels.end() ){
            --itA;
          }

        }
        else{
          ++itB;
          while(itB != bLevels.end() ){
            rB = getMeanRadius(*itB, p);
            TPrecision tmp = rA-rB;
            if( fabs(tmp) > fabs(delta) ){
              --itB;
              break;
            }
            else{
              delta = tmp;
            }
            ++itB;
          } 

          if(itB == bLevels.end() ){
            --itB;
          }
        }

        
      }


      MultiscaleTransportLevel<TPrecision> *A = *itA;
      MultiscaleTransportLevel<TPrecision> *B = *itB;
      
      TransportPlanSolutions<TPrecision> *prevSol = NULL;
      int currentScale = 0;
      while( itA != aLevels.end() || itB != bLevels.end() ){

	std::cout << std::endl << std::endl << "Scale: " << currentScale;
        std::cout << std::endl << std::endl;
        ++currentScale;

	//Check if this is the last scale
        bool lastScale = false;
        if(itA == aLevels.end() ){
          typename std::vector< MultiscaleTransportLevel<TPrecision> * >::iterator
            it = itB;
          it++;
          lastScale = it == bLevels.end();
        }
        else if(itB == bLevels.end() ){
          typename std::vector< MultiscaleTransportLevel<TPrecision> * >::iterator
            it = itA;
          it++;
          lastScale = it == aLevels.end();
        }
        else{
         typename std::vector< MultiscaleTransportLevel<TPrecision> * >::iterator
            it = itA;
          it++;
          lastScale = it == aLevels.end();
          it = itB;
          it++;
          lastScale = lastScale && it == bLevels.end();
        }


        TransportPlanSolutions<TPrecision> *sol = solveLP(A, B, prevSol, p, lastScale);
        solutions.push_back( sol->getPrimarySolution() );
        prevSol = sol;
       

        if(itA != aLevels.end()){
         ++itA;
        } 
        if( itA != aLevels.end()){
          A = *itA;
        }

        if( itB != bLevels.end()){
          ++itB;
        }
        if( itB != bLevels.end()){
          B = *itB;
        }

      }



      return solutions;
    };






  protected:


    virtual TransportPlanSolutions<TPrecision>
      *solveLP(MultiscaleTransportLevel<TPrecision> *source,
          MultiscaleTransportLevel<TPrecision> *target, TransportPlanSolutions<TPrecision>
          *prevSol, double p, bool lastScale) = 0;


    MultiscaleTransportLevel<TPrecision> *getRootSource(){
      return rootSource;
    };

    MultiscaleTransportLevel<TPrecision> *getRootTarget(){
      return rootTarget;
    };



  private:

    MultiscaleTransportLevel<TPrecision> *rootSource;
    MultiscaleTransportLevel<TPrecision> *rootTarget;

    void normalize(std::vector< MultiscaleTransportLevel<TPrecision> * > &levels){
      //Normalize masses to one at each level and match up parent child mass
      //relations
      for(typename std::vector< MultiscaleTransportLevel<TPrecision> *
          >::reverse_iterator it =
        levels.rbegin(); it != levels.rend(); ++it){
        TransportNodeVector &n = (*it)->getNodes();
        TPrecision total = 0;
        for(TransportNodeVectorCIterator nIt = n.begin(); nIt != n.end(); nIt++){
          const TransportNodeVector &kids = (*nIt)->getChildren();
          TPrecision sum = 0;
          //if(kids.size() == 0 ){
          //  sum = (*nIt)->getMass();
          //}
          for(TransportNodeVectorCIterator kIt = kids.begin(); kIt !=kids.end();
              ++kIt){
            sum += (*kIt)->getMass();
          }
          //std::cout << (*nIt)->getMass() - sum << std::endl;
          (*nIt)->setMass(sum);
          total += sum;
        }
        for(TransportNodeVectorIterator nIt = n.begin(); nIt != n.end(); nIt++){
          (*nIt)->setMass((*nIt)->getMass()/total);
        }
      }
    };



    TPrecision getMeanRadius(MultiscaleTransportLevel<TPrecision> *level, double p){
        TransportNodeVector &n = level->getNodes();
        TPrecision radius = 0;
        for(TransportNodeVectorCIterator nIt = n.begin(); nIt != n.end(); nIt++){
          radius += (*nIt)->getLocalNodeRadius();
        }
        return radius / n.size();
    };

   
};


#endif
