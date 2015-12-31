#ifndef IPCANODEFACTORY_H
#define IPCANODEFACTORY_H


#include "IPCANode.h"
#include "EigenRandomSVD.h"

#include <queue>
#include <map>

#define STREAMING_THRESHOLD 100000000



//Node factory decides on the dimensionality of each node
template <typename TPrecision>
class IPCANodeFactory{

  public:



    typedef typename GMRANode<TPrecision>::MatrixXp MatrixXp;
    typedef typename GMRANode<TPrecision>::VectorXp VectorXp;

    enum StoppingCriterium {TOTAL_R2, NODE_R2, NODE_MSE, NODE_RADIUS,
      RELATIVE_NODE_RADIUS, MASS_RADIUS};

    enum SplitDirectionStrategy {PC, RANDOM_PC, AXIS_ALIGNED, LDA};
    enum SplitStrategy{MEAN, MIDPOINT, MEDIAN, RANDOM_MEAN, RANDOM_MIDPOINT};
    
    StoppingCriterium stop;
    TPrecision epsilon;
    SplitStrategy splitStrategy;
    unsigned int minPoints;
    SplitDirectionStrategy splitDirectionStrategy;
    int maxKidDim;

  protected:

    int maxDim;
    TPrecision radius;
    TPrecision totalVar;
    VectorXp sigma;
    MatrixXp phi;
    GMRADataObject<TPrecision> *data;
    
    virtual void truncateSVD() = 0;

  private:
    


    void computeSVD(std::vector<int> &indices, const VectorXp &mean){      

      using namespace Eigen;
      
      //Sample randomly from if matrix too large  
      MatrixXp X;
      int d = data->dimension();
      if( (indices.size() * d) > STREAMING_THRESHOLD){
        int n = STREAMING_THRESHOLD / d;
        X = MatrixXp( d, n );
        for(unsigned int i=0; i < X.cols(); i++){
          static Random<TPrecision> random;
          int index = std::floor(random.Uniform() * indices.size() );
          X.col(i)  = data->getPoint( indices[index] ) - mean;
        }

      }
      else{
        X = MatrixXp( d, indices.size() );        
        for(unsigned int i=0; i < X.cols(); i++){
          X.col(i)  = data->getPoint( indices[i] ) - mean;
        }
      }


      totalVar = 0;
      radius = 0;
      for(unsigned int i=0; i < X.cols(); i++){
        TPrecision tmp = X.col(i).squaredNorm();
        totalVar += tmp;
        radius = std::max(radius, tmp );
      }
      if(radius > 0 ){
        radius = sqrt(radius);
      }
      if( X.cols() > 1){
        totalVar /= (X.cols()-1);
      }


      //Select method for computing svd depending on size of matrix

      //trivial case
      if(X.cols() == 1){
        phi = MatrixXp::Zero(this->data->dimension(), 1);
        sigma = VectorXp::Zero(1);
      }
      //Do normal svd for smaller dimensions
      else if(X.cols() < 1000 && X.rows() < maxDim+20  ){        
        JacobiSVD<MatrixXd> svd(X, ComputeThinU | ComputeThinV);
        sigma = svd.singularValues();
        phi = svd.matrixU();
      }
      //else do randomized svd
      else{ 
        EigenLinalg::RandomSVD<TPrecision> svd(X, maxDim+4, 1);
        phi = svd.GetU();
        sigma = svd.GetS();
      }
      if(X.cols() > 1){
        sigma.array() /= sqrt( X.cols()-1.0 );
      }
      
      this->truncateSVD();

    };






  public:
    
    
    IPCANodeFactory(GMRADataObject<TPrecision> *D, int maxD):data(D){       
      
      maxDim        = std::max(1, maxD); 
      stop          = TOTAL_R2;
      epsilon       = 0;
      splitStrategy = MEAN;
      minPoints     = 1;
      splitDirectionStrategy = PC;
      maxKidDim     = 3;

    };
    
    
    virtual ~IPCANodeFactory(){
    };


    GMRANode<TPrecision> *createNode(std::vector<int> &indices){

      VectorXp mean = VectorXp::Zero( this->data->dimension() );
      for( int i=0; i< indices.size(); i++ ){
        mean += this->data->getPoint( indices[i] );
      }
      mean.array() /= indices.size();

      computeSVD(indices, mean);


      MatrixXp dir;
      //Split direction
      if(splitDirectionStrategy == PC){
        dir =phi.leftCols( std::min( maxKidDim,  (int) phi.cols() )  );
      }
      else if(splitDirectionStrategy == RANDOM_PC){
        static Random<TPrecision> random;
        dir = MatrixXp::Zero( phi.rows(), std::min( maxKidDim, (int) phi.cols() )  );
        for(int k=0; k<dir.cols(); k++){
          for(unsigned int i=0; i<dir.cols(); i++){
            double w = random.Uniform();
            dir.col(k) += w*phi.col(i);
          }
          dir.col(k).normalize();
        }
      }
      else if(splitDirectionStrategy == AXIS_ALIGNED){
        dir = MatrixXp::Zero(phi.rows(), phi.cols());
        for(int i=0;i<dir.cols(); i++){
          dir(i,i) = 1;
        }
      }



      //Split location
      VectorXp splitCenter;
      if(splitStrategy == MEAN || 
         splitStrategy == RANDOM_MEAN){

        splitCenter = mean;
        if(splitStrategy == RANDOM_MEAN){
          for(int i=0; i < dir.cols() ; i++){
            static Random<TPrecision> random;
            TPrecision s = random.Normal() * sigma(i);
            splitCenter += dir.col(i) * s;
          }
        }
      }
      else if(splitStrategy == MIDPOINT || 
              splitStrategy == RANDOM_MIDPOINT){
       //TODO: Broken, wrong, bad, etc. 
        VectorXp minP = -5 * sigma;
        VectorXp maxP = 5 * sigma;
        VectorXp mid = minP;
        
        if(splitStrategy == RANDOM_MIDPOINT){
          static Random<TPrecision> random;
          for(int i=0; i< mid.size(); i++){
            TPrecision s = random.Uniform()*0.4 + 0.3;
            mid(i) += ( maxP(i)-minP(i) * s );
          }
        }
        else{         
          mid += maxP;
          mid *= 0.5; 
        }
        splitCenter = dir * mid + mean;

      }
      else if(splitStrategy == MEDIAN){

      }


      VectorXp a = dir.transpose() * splitCenter;


      IPCANode<TPrecision> *node = new IPCANode<TPrecision>(mean, indices, phi, sigma, radius, dir, a, totalVar);

      return node;
    };
  


    bool isStop(IPCANode<TPrecision> *node, int nPoints, TPrecision rootRadius, TPrecision rootMSE, TPrecision rootVariance){

      VectorXp sigma2 = node->getSigma().array().square(); 
      TPrecision mse = node->getTotalVariance() - sigma2.sum();
      TPrecision R2 = 1 - mse / node->getTotalVariance() ;
      TPrecision relativeR2 = 1- mse / rootVariance;
      TPrecision relativeRadius = node->getL2Radius() / rootRadius;
      TPrecision mass = node->getPoints().size() / (TPrecision) nPoints;
      TPrecision massRadius = mass * node->getL2Radius();



#ifdef VERBOSE
      std::cout << "Node Variance : " << node->getTotalVariance() << std::endl;
      std::cout << "Node MSE : " << mse << std::endl;
      std::cout << "Node sigma : " << node->getSigma() << std::endl;
      std::cout << "MinPoints : " << minPoints << std::endl;
      std::cout << "Node size : " << node->getPoints().size() << std::endl;
      std::cout << "Node R^2 : " << R2 << std::endl;
      std::cout << "Node total R^2 : " << relativeR2 << std::endl;
      std::cout << "Node L2 radius : " << node->getL2Radius() << std::endl;
      std::cout << "Node relative L2 radius : " << relativeRadius << std::endl;
      std::cout << "Node mass : " << mass << std::endl;
      std::cout << "Node mass * radius : " << massRadius << std::endl;
      std::cout << "Node dim : " << node->getPhi().cols() << std::endl << std::endl;
#endif


      if( node->getPoints().size() <= std::max(minPoints, (unsigned int) 1 ) ){ 
        return true;
      }
      if(TOTAL_R2 && relativeR2 >= epsilon){
        return true;
      }
      if(stop == NODE_R2 && R2 >= epsilon){
        return true;
      }
      if(stop == NODE_MSE && mse <= epsilon){
        return true;
      }
      if(stop == NODE_RADIUS && node->getL2Radius() <= epsilon){
        return true;
      }
      if(stop == RELATIVE_NODE_RADIUS && relativeRadius <= epsilon){
        return true;
      }      
      if(stop == MASS_RADIUS && massRadius <= epsilon){
        return true;
      }

      return false;
    };


    void setSplitStrategy(int s){
      switch(s){
        case 1:
          std::cout << "Mean" << std::endl;
          splitStrategy = MEAN;
        case 2:
          std::cout << "Mid" << std::endl;
          splitStrategy = MIDPOINT;
        case 3:
          std::cout << "Random Mean" << std::endl;
          splitStrategy =  RANDOM_MEAN;
        case 4:
          std::cout << "Random Mid" << std::endl;
          splitStrategy = RANDOM_MIDPOINT;
      }

      splitStrategy = MEAN;

    };

    void setStoppingCriterium(int s){
      switch(s){
        case 1:
          std::cout << "Total R^2" << std::endl;
          stop = TOTAL_R2;
        case 2:
          std::cout << "Node R^2" << std::endl;
          stop = NODE_R2;
        case 3:
          std::cout << "Node MSE" << std::endl;
          stop = NODE_MSE;
        case 4:
          std::cout << "Node Radius" << std::endl;
          stop = NODE_RADIUS;
        case 5:
          std::cout << "Relative Node Radius" << std::endl;
          stop = RELATIVE_NODE_RADIUS;      
        case 6:
          std::cout << "Mass Radius" << std::endl;
          stop = MASS_RADIUS;      
      }
      stop = NODE_RADIUS;
    };



   void setSplitDirectionStrategy(int s){
      switch(s){
        case 1:
          splitDirectionStrategy = PC;
        case 2:
          splitDirectionStrategy = RANDOM_PC;
        case 3:
          splitDirectionStrategy = AXIS_ALIGNED;
      }
      splitDirectionStrategy= PC;
    };


};








template <typename TPrecision>
class FixedNodeFactory : public IPCANodeFactory<TPrecision>{
   public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;

  protected:
   void truncateSVD(){     
     int d = this->maxDim; 
     if(this->phi.cols() > d){
        MatrixXp phiTmp = this->phi.leftCols(d);
        this->phi = phiTmp;
        VectorXp tmp = this->sigma.head(d);
        this->sigma = tmp;
     }
   };


  public:
    FixedNodeFactory(GMRADataObject<TPrecision> *D, int dim = 5) : IPCANodeFactory<TPrecision>(D, dim){
    };

};





template <typename TPrecision>
class RelativePrecisionNodeFactory : public IPCANodeFactory<TPrecision>{
  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
  private:
    TPrecision t;  
  
  protected:
    void truncateSVD(){        
      
      int d = 1;
      VectorXp sigma2 = this->sigma.array().square();
      TPrecision mse0 = this->totalVar;
      TPrecision mseTmp = 0;
      for(int i=0; i<sigma2.size(); i++){
        mseTmp += sigma2(i);
        if(mseTmp/mse0 > t){
          break;
        }
        d++;
        if(d == this->maxDim){
          break;
        }
      }

   
       if( d < this->phi.cols() ){
        MatrixXp phiTmp = this->phi.leftCols(d);
        this->phi = phiTmp;
        VectorXp tmp = this->sigma.head(d);
        this->sigma = tmp;
      };
    };


  public:
    RelativePrecisionNodeFactory(GMRADataObject<TPrecision> *D, int maxDim, TPrecision rel) : IPCANodeFactory<TPrecision>(D, maxDim), t(rel){
    };
    


};


template <typename TPrecision>
class RelativeRatioNodeFactory : public IPCANodeFactory<TPrecision>{
  public:
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
  private:
    TPrecision t;  
  
  
  protected:
    void truncateSVD(){        

      int d = 1;
      VectorXp sigma2 = this->sigma.array().square();
      for(int i=1; i<this->sigma.size(); i++){
        if(sigma2(i-1)/sigma2(i) < t){
          break;
        }
        d++;
        if(d == this->maxDim){
          break;
        }
      }
   
      if(this->phi.cols() > d){
        MatrixXp phiTmp = this->phi.leftCols(d);
        this->phi = phiTmp;
        VectorXp tmp = this->sigma.head(d);
        this->sigma = tmp;
      }

    };


  public:
    RelativeRatioNodeFactory(GMRADataObject<TPrecision> *D, int maxDim, TPrecision ratio) :
      IPCANodeFactory<TPrecision>(D, maxDim), t(ratio){ 
        
    };

   
};









#endif
