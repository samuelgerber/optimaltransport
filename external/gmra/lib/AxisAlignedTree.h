#ifndef AATREE_H
#define AATREE_H

#include "GMRATree.h"
#include "SVD.h"
#include "Random.h"


#include <map>
#include <iostream>
#include <fstream>





//GMRANode subclass
template <typename TPrecision> 
class AANode : public GMRANodeBase<TPrecision>{
  public:
    typedef typename GMRANode<TPrecision>::NodeVector NodeVector;
    typedef typename NodeVector::iterator NodeVectorIterator;


  private:
    NodeVector children;
    std::map<int, int> childmap;
    TPrecision mass;

  public:
    enum SplitStrategy {MEAN, MIDPOINT, MEDIAN, RANDOM_MEAN, RANDOM_MIDPOINT};

    FortranLinalg::DenseMatrix<TPrecision> dir;
    FortranLinalg::DenseVector<TPrecision> a; 



    virtual ~AANode(){
      dir.deallocate();    
      center.deallocate();
      a.deallocate();
    };



    AANode<TPrecision>(){
    };



    AANode<TPrecision>(FortranLinalg::DenseMatrix<TPrecision> &Xi,
        FortranLinalg::DenseVector<TPrecision> w, std::vector<int> &ind,
        SplitStrategy splitPoint = MEAN, int dim = 3) :
      indices(ind){

      using namespace FortranLinalg;


      DenseMatrix<TPrecision> Xc = createMatrix(Xi);
      
      center = DenseVector<TPrecision>( Xc.M() ) ;
      Linalg<TPrecision>::Zero(center);
      mass = 0;
      for(int i=0; i<ind.size(); i++){
        int index = ind[i];

        for(int k= 0; k<center.N(); k++){
          center(k) += Xi(k, index)*w(index);
        }
        mass += w(index); 
      }

      Linalg<TPrecision>::Scale(center, 1.0/mass, center);
      Linalg<TPrecision>::SubtractColumnwise(Xc, center, Xc);
      
      DenseVector<TPrecision> l = Linalg<TPrecision>::ColumnwiseSquaredNorm(Xc);
      l2Radius =  sqrtf(Linalg<TPrecision>::Max(l));
      l.deallocate();

      dir = DenseMatrix<TPrecision>(Xc.M(), dim);
      Linalg<TPrecision>::Zero(dir);
      for(int i=0;i<dim; i++){
        dir(i,i) = 1;
      }

      DenseVector<TPrecision> splitCenter;
      if(splitPoint == MEAN || splitPoint == RANDOM_MEAN){
        splitCenter = Linalg<TPrecision>::Copy(center);
        if(splitPoint == RANDOM_MEAN){
          for(int i=0; i < dir.N() ; i++){
            DenseVector<TPrecision> tmp = Linalg<TPrecision>::ExtractColumn(dir, i);
            static Random<TPrecision> random;
            TPrecision s = random.Normal() * l2Radius/3.0;
            Linalg<TPrecision>::Scale(tmp, s, tmp);
            Linalg<TPrecision>::Add(splitCenter, tmp, splitCenter);
            tmp.deallocate();
          }
        }
      }
      else if(splitPoint == MIDPOINT || splitPoint ==  RANDOM_MIDPOINT){
        
        DenseMatrix<TPrecision> Xp = Linalg<TPrecision>::Multiply(dir, Xc, true);
        DenseVector<TPrecision> minP = Linalg<TPrecision>::RowMin(Xp);
        DenseVector<TPrecision> maxP = Linalg<TPrecision>::RowMax(Xp);
        DenseVector<TPrecision> mid = Linalg<TPrecision>::Copy(minP);
        
        if(splitPoint == RANDOM_MIDPOINT){
          static Random<TPrecision> random;
          for(int i=0; i< mid.N(); i++){
            TPrecision s = random.Uniform()*0.4 + 0.3;
            mid(i) += ( maxP(i)-minP(i) * s );
          }
        }
        else{          
          Linalg<TPrecision>::Add(mid, maxP, mid);
          Linalg<TPrecision>::Scale(mid, 0.5, mid);
        }
        splitCenter = Linalg<TPrecision>::Multiply(dir, mid);
        Linalg<TPrecision>::Add(splitCenter, center, splitCenter);

        Xp.deallocate();
        mid.deallocate();
        maxP.deallocate();
        minP.deallocate();

      }
      else if(splitPoint == MEDIAN){
        /*
           DenseVector<TPrecision> a = Linalg<TPrecision>::Min(Vt);
           DenseVector<TPrecision> b = Linalg<TPrecision>::Max(Vt);
           Linalg<TPrecision>::Subtract(b, a, b);
           Linalg<TPrecision>::Scale(b, 0.5, b);
           Linalg<TPrecision>::Add(a, b, a);
           splitCenter = Linalg<TPrecision>::Multiply(U, a);
           Linalg<TPrecision>::Add(splitCenter, center, splitCenter);
           a.deallocate();
           b.deallocate();
         */
      }


      a = Linalg<TPrecision>::Multiply(dir, splitCenter, true);
      splitCenter.deallocate();

      Xc.deallocate();



    };


    TPrecision getDimensionFreeVolume(){
      return this->getRadius(); 
    };



    void addChild(GMRANode<TPrecision> *node, int childIndex){
      children.push_back(node);
      childmap[childIndex] = children.size()-1;
    };




    GMRANode<TPrecision> *getChild(int childIndex){
      std::map<int, int>::iterator it = childmap.find(childIndex);
      if(it == childmap.end()){
        return NULL;
      }
      return children[it->second];
    };




    int getChildIndex(FortranLinalg::DenseVector<TPrecision> x){
      using namespace FortranLinalg;
      DenseVector<TPrecision> s = Linalg<TPrecision>::Multiply(dir, x, true);
      int childIndex = 0;
      int factor=1;
      for(int i=0; i<s.N(); i++){
        if( s(i) > a(i) ){
          childIndex += factor;
        }
        factor *= 2;
      }
      s.deallocate();
      return childIndex;

    };



    virtual GMRANode<TPrecision> *findDescendant( FortranLinalg::DenseVector<TPrecision> &x ){
       int index = getChildIndex(x);
       return getChild(index);
    };




    virtual NodeVector &getChildren(){
       return children;
    };


    std::vector<int> &getPoints(){
      return indices;
    };

    FortranLinalg::DenseVector<TPrecision> &getCenter(){
      return center;
    };


    TPrecision getMass(){
      return mass;
    };

    TPrecision getL2Radius(){
      return l2Radius;
    };

    void flatten(std::ostream &file){
      int nSplit = dir.N();
      file.write((char*)&nSplit, sizeof(int) );
      
      int nKids = children.size();
      file.write((char*)&nKids, sizeof(int) );

      for(std::map<int, int>::iterator it = childmap.begin(); it !=
          childmap.end(); it++){
        file.write((char*)&it->first , sizeof(int) );
        file.write((char*)&it->second, sizeof(int) );
      }

      file.write((char*)center.data(), center.N() * sizeof(TPrecision) );
      file.write((char*)dir.data(), dir.N() * dir.M()*sizeof(TPrecision) );
      file.write((char*)a.data(), a.N() * sizeof(TPrecision) );

      int nPoints = indices.size();
      file.write((char*) &nPoints, sizeof(int) );
      file.write((char*) indices.data(), indices.size()*sizeof(int) );

      file.write((char*) &l2Radius, sizeof(TPrecision) );
      file.write((char*) &mass, sizeof(TPrecision) );

    };




    int unflatten(std::istream &file, int m){
      using namespace FortranLinalg;
      
      int nSplit = 0;
      file.read( (char*) &nSplit, sizeof(int) );
      
      int nKids = 0;
      file.read( (char*) &nKids, sizeof(int) );

      for(int i=0; i<nKids; i++){
        int i1 = 0;
        file.read( (char*) &i1, sizeof(int) );
        int i2 = 0;
        file.read( (char*) &i2, sizeof(int) );
        childmap[i1] = i2;
      }

      

      center = DenseVector<TPrecision>(m);
      file.read((char*)center.data(), m*sizeof(TPrecision));

      dir = DenseMatrix<TPrecision>(m, nSplit);
      file.read((char*)dir.data(), m*nSplit*sizeof(TPrecision));

      a = DenseVector<TPrecision>(nSplit);
      file.read((char*)a.data(), sizeof(TPrecision)*nSplit );

      int nPoints;
      file.read( (char*) &nPoints, sizeof(int) );
      indices.resize(nPoints);
      file.read((char*)indices.data(), nPoints*sizeof(int) );

      file.read((char*) &l2Radius, sizeof(TPrecision) );
      file.read((char*) &mass, sizeof(TPrecision) );



      return nKids;
    };




    virtual void translate(FortranLinalg::DenseVector<TPrecision> &x){
      using namespace FortranLinalg;
      Linalg<TPrecision>::Add( center, x, center );
    };




    virtual void affine(FortranLinalg::DenseMatrix<TPrecision> &A){
      using namespace FortranLinalg;
     
      DenseVector<TPrecision> centerN = Linalg<TPrecision>::Multiply(A, center);
      center.deallocate();
      center = centerN;
      
    };



    FortranLinalg::DenseVector<TPrecision> center;
    TPrecision l2Radius;


    std::vector<int> indices;



  private:

    FortranLinalg::DenseMatrix<TPrecision> createMatrix(FortranLinalg::DenseMatrix<TPrecision> &X){
      FortranLinalg::DenseMatrix<TPrecision> M(X.M(), indices.size());

      std::vector<int>::iterator it = indices.begin();
      for(unsigned int i=0; i < M.N(); i++, ++it){
        FortranLinalg::Linalg<TPrecision>::SetColumn(M, i, X, *it);
      }

      return M;
    };

};






template <typename TPrecision>
class AxisAlignedTree : public GMRATree<TPrecision>{

  public:
    enum StoppingCriterium {NODE_RADIUS, RELATIVE_NODE_RADIUS, MASS_RADIUS};

  private:
    typedef typename AANode<TPrecision>::SplitStrategy SplitStrategy;

    StoppingCriterium stop;
    TPrecision epsilon;
    SplitStrategy splitStrategy;
    unsigned int minPoints;
    FortranLinalg::DenseMatrix<TPrecision> Xref;
    int dim;
    FortranLinalg::DenseVector<TPrecision> weights;





    void buildTreeRecursive(AANode<TPrecision> *node, FortranLinalg::DenseMatrix<TPrecision>
        &X, FortranLinalg::DenseVector<TPrecision> &w, TPrecision rootRadius, TPrecision nPoints){

      using namespace FortranLinalg;




      TPrecision relativeRadius = node->getL2Radius() / rootRadius;
      TPrecision mass = node->getPoints().size() / nPoints;
      TPrecision massRadius = mass * node->getL2Radius();



#ifdef VERBOSE
      std::cout << "MinPoints : " << minPoints << std::endl;
      std::cout << "Node L2 radius : " << node->getL2Radius() << std::endl;
      std::cout << "Node relative L2 radius : " << relativeRadius << std::endl;
      std::cout << "Node mass : " << mass << std::endl;
      std::cout << "Node mass * radius : " << massRadius << std::endl;
#endif


      if( node->getPoints().size() <= std::max(minPoints, (unsigned int) 1 ) ){ 
        return;
      }
      if(stop == NODE_RADIUS && node->getL2Radius() <= epsilon){
        return;
      }
      if(stop == RELATIVE_NODE_RADIUS && relativeRadius <= epsilon){
        return;
      }      
      if(stop == MASS_RADIUS && massRadius <= epsilon){
        return;
      }

      int size = pow( (double) 2, dim);
      std::vector< std::vector<int> > children(size);
      DenseVector<TPrecision> tmp(X.M());
      std::vector<int> &nodePts = node->getPoints();
      for(std::vector<int>::iterator it = nodePts.begin();
          it!=nodePts.end(); ++it){
        int i = *it;
        Linalg<TPrecision>::ExtractColumn(X, i, tmp);
        int childIndex = node->getChildIndex(tmp);
        children[childIndex].push_back(i);
      }
      tmp.deallocate();

      for(int i=0; i< children.size(); i++){
        if(children[i].size() > 0){
          AANode<TPrecision> *n = new AANode<TPrecision>(X, w, children[i],
              splitStrategy, dim);
          node->addChild(n, i);
          buildTreeRecursive(n, X, w, rootRadius, nPoints);
        }
      }

    };





  public:


    AxisAlignedTree(){};



    AxisAlignedTree(TPrecision eps, StoppingCriterium stopC,  SplitStrategy ss,
        int d = 3, int minLeafSize=1): stop(stopC), epsilon(eps),
    splitStrategy(ss), minPoints(minLeafSize), dim(d){

    };
    
    AxisAlignedTree(FortranLinalg::DenseVector<TPrecision> w, TPrecision eps,
        StoppingCriterium stopC,  SplitStrategy ss, int d = 3, int
        minLeafSize=1): stop(stopC), epsilon(eps), splitStrategy(ss),
    minPoints(minLeafSize), dim(d), weights(w){

    };




    virtual ~AxisAlignedTree(){
      delete[] Xref.getColumnAccessor();
    };




    //
    void construct(TPrecision *X, int m, int n){

      FortranLinalg::DenseMatrix<TPrecision> Xd(m, n, X);
      Xref = Xd;

      std::vector<int> all;
      for(int i=0; i<n; i++){
        all.push_back(i);
      }

      FortranLinalg::DenseVector<TPrecision> w;
      if(weights.N() == 0){
        w = FortranLinalg::DenseVector<TPrecision>(n);
        FortranLinalg::Linalg<TPrecision>::Set(w, 1);
      }
      else{
        w = weights;
      }

      
      AANode<TPrecision> *root = new AANode<TPrecision>(Xd, w, all, splitStrategy, dim);

      buildTreeRecursive( (AANode<TPrecision>*) root, Xd, w, root->getL2Radius(),
          root->getPoints().size()); 

      this->setRoot(root);
      this->setupParents();
    
      if(weights.N() == 0){
        w.deallocate();
      }

    };




    FortranLinalg::DenseVector<TPrecision> getPoint(int index){
      FortranLinalg::DenseVector<TPrecision> x = FortranLinalg::Linalg<TPrecision>::ExtractColumn(Xref, index);
      return x;
    };




    void add(TPrecision *x){
      throw "not implemented yet";
    };




    std::vector<GMRANode<TPrecision> *> getLeafPath( double *x ) {
      using namespace FortranLinalg;
   
      GMRANode<TPrecision> *node = this->getRoot();
      DenseVector<TPrecision> xv( node->getCenter().N(), x);

      std::vector<GMRANode<TPrecision> *> path;
      while(node != NULL ){
        path.push_back( node );
        node = node->findDescendant( xv );
      }

      return path;
    }; 






    void flatten( std::ofstream &file ){
      using namespace FortranLinalg;
  
      std::list<AANode<TPrecision> *> nodes;
      AANode<TPrecision>* root = (AANode<TPrecision>*) this->getRoot();
      nodes.push_back(root);

      file.write( (char*) &epsilon, sizeof(TPrecision) ) ;
      unsigned int m = dim;
      file.write( (char*) &m, sizeof(unsigned int) ) ;   
      file.write( (char*) &minPoints, sizeof(unsigned int) );   

      while( !nodes.empty() ){
        AANode<TPrecision> *node =  nodes.front();
        nodes.pop_front();

        node->flatten(file);

        for(typename std::vector< GMRANode<TPrecision> * >::iterator it =
            node->getChildren().begin(); it != node->getChildren().end(); ++it){
          nodes.push_back((AANode<TPrecision>*) *it);
        }
      }

    };






    void unflatten( std::ifstream &file ){
      using namespace FortranLinalg;
    
      file.read( (char*) &epsilon, sizeof(TPrecision) ) ;
      unsigned int m;   
      file.read( (char*) &m, sizeof(unsigned int) );   
      dim = m;

      file.read( (char*) &minPoints, sizeof(unsigned int) );   

      std::list<AANode<TPrecision> *> nodes;
      std::list<int> nKids;

      AANode<TPrecision> *cur = NULL;
      int nK = 0;
      while( !file.eof() ){
        
        AANode<TPrecision> *node = new AANode<TPrecision>();
        int n = node->unflatten(file, m);
        if(cur == NULL){
          this->setRoot(node);
          cur = node;
          nK = n;
        }
        else{ 
          if( n > 0 ){
            nodes.push_back(node);
            nKids.push_back(n);
          }
          if( cur->getChildren().size() == nK ){
            cur = nodes.front();
            nodes.pop_front();

            nK = nKids.front();
            nKids.pop_front();
          }
          cur->getChildren().push_back(node);

        }
        file.peek();
      }

      this->setupParents();

    };


};


#endif
