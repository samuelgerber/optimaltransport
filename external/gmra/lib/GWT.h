#ifndef GWT_H
#define GWT_H


#include "EigenRandomSVD.h"

#include <algorithm>
#include <map>

#include "GMRATree.h"
//Interface for GWT

template <typename TPrecision> 
class GWTNode : public GMRANodeDecorator<TPrecision>{
  public:  
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
  
  private:
    VectorXp w;

    MatrixXp psi;
    MatrixXp psi_t_phi;
   
    int ID;
  

  public:


    //Caller is responsible for cleaning up passed in variables properly
    //This class stores pointers to cx, s, T and p
    GWTNode<TPrecision>(const MatrixXp &phi, const VectorXp &center, int id,
        GWTNode<TPrecision> *p, GMRANode<TPrecision> *node ) :
      GMRANodeDecorator<TPrecision>(node) {
        
      
      ID = id;
      //phi_t_phi = Linalg<TPrecision>::Multiply(phi, phi, true);

      if(p != NULL){
        MatrixXp &p_phi = p->getPhi();

        MatrixXp phi_p = p_phi.transpose() * phi;
        psi = phi - p_phi * phi_p;
        
        Eigen::JacobiSVD<MatrixXp> svd(psi, Eigen::ComputeThinU | Eigen::ComputeThinV);

        int cut = svd.singularValues().size();
        for(int i=0; i<svd.singularValues().size(); i++){
          if( svd.singularValues()(i) < 0.00000000001 ){
            cut = i;
            break;
          }
        }
        if(cut == 0){
          psi = MatrixXp(0, 0);
        }
        else{
          psi = svd.matrixU().leftCols(cut);
        }

        
        if(cut == 0){
          psi_t_phi = MatrixXp(0, 0);
        }
        else{
          psi_t_phi = psi.transpose() * phi;
        }
        VectorXp t = center - p->getCenter();
        w = t - p_phi * (p_phi.transpose() * t );

      }

    };


    virtual  ~GWTNode(){
    };



    MatrixXp &getPsi(){
      return psi;
    };

    MatrixXp &getPsitPhi(){
      return psi_t_phi;
    };


    VectorXp &getW(){
      return w;
    };

    int getID(){
      return ID;
    };

    virtual MatrixXp &getPhi() = 0;
    virtual VectorXp &getSigma() = 0;
    virtual VectorXp &getCenter() = 0;

};






template <typename TPrecision> 
class GWTCoefficients{

  public:  
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
    
    GWTCoefficients(){
      maxD = 0;
    };

    //Matrix of coefficents at each 
    //scale from finest (0) to coarsest (coeff.M()-1)
    std::vector< VectorXp > coeff;
    std::vector< int > ids;
    int maxD;

    //Node ids for each scale
    GMRANode<TPrecision> *leaf;
};





template <typename TPrecision>
class GWT{

  public:  
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;
  
  
  private:
  
    GMRATree<TPrecision> *tree;


  protected:
    //Subclass needs to decorate ech GMRATree node info to each GMRANode in the tree.
    virtual void setupGWT(GMRATree<TPrecision> *tree) = 0;

    
  public:

    
    //Setup tree for GWT, i.e. add covariance structur to each node -> see
    //setupGWT
    void setTree(GMRATree<TPrecision> *gTree){
      tree = gTree;
      setupGWT(tree);
    };

    

    //Encoding of GMRA coefficents 
    GWTCoefficients<TPrecision> encode(const VectorXp &x){
      GMRANode<TPrecision> *node = tree->getLeafPath( x ).back();
      GWTNode<TPrecision> *gwt = dynamic_cast<GWTNode<TPrecision> *>( node );
      GWTNode<TPrecision> *parent = dynamic_cast<GWTNode<TPrecision> *>(
          node->getParent() );
      
      GWTCoefficients<TPrecision> c;
      c.leaf = node;

      std::list< VectorXp > ctmp;


      VectorXp delta = x - gwt->getCenter();
      VectorXp p = gwt->getPhi().transpose() * delta;
      VectorXp x_J = gwt->getPhi() * p;
      x_J += gwt->getCenter();

      while(parent != NULL){
        VectorXp q;
        if(gwt->getPsi().cols() != 0){
          q = gwt->getPsitPhi() * p;
        }
        c.coeff.push_back( q );
        c.ids.push_back( gwt->getID() );
        if(q.N() > c.maxD){
          c.maxD = q.N();
        }

        delta = x_J - parent->getCenter();
        p = parent->getPhi().transpose() * delta;

        gwt = parent;
        parent = dynamic_cast< GWTNode<TPrecision> *>( gwt->getParent() );
      }

      c.coeff.push_back(p);
      c.ids.push_back(gwt->getID());
      std::reverse(c.coeff.begin(), c.coeff.end() );
      std::reverse(c.ids.begin(), c.ids.end() );
      if(p.N() > c.maxD){
        c.maxD = p.N();
      }

      return c;
    };




    //Decoding from GMRA coefficents at specific scale (0 = coarsest,
    //a.coeff.M() =finest also -1 = finest)
    VectorXp decode(GWTCoefficients<TPrecision> &a, int scale =-1){
      GMRANode<TPrecision> *gmraNode = a.leaf; 
      GWTNode<TPrecision> *gwt = (GWTNode<TPrecision> *) gmraNode;
      GWTNode<TPrecision> *parent = (GWTNode<TPrecision> *)
        gmraNode->getParent();
      

      int nScales = 0;
      GWTNode<TPrecision> *gwtTmp = gwt;
      while(gwtTmp->getParent() != NULL){
        nScales++;
        gwtTmp = dynamic_cast< GWTNode<TPrecision> *>( gwtTmp->getParent() );
      }
       
      if(scale < 0 || scale > nScales){
        scale = nScales;
      }
      int scaleIndex = nScales;
      while(scaleIndex != scale && gwt->getParent() != NULL){
          scaleIndex--;
          gwt = dynamic_cast< GWTNode<TPrecision> *>( gwt->getParent() );
      }
        
      
      VectorXp sumq;  
      if(parent != NULL){
        if(gwt->getPsi().cols() != 0){
          VectorXp q = a.coeff[scaleIndex];
          sumq = gwt->getPsi() * q.head(gwt->getPsi().cols()); 
          sumq += gwt->getW(); 
        }
        else{
          sumq = gwt->getW(); 
        }

        
        gwt = parent;
        parent = (GWTNode<TPrecision> *) gwt->getParent();

        scaleIndex--;
        while(parent != NULL){
          VectorXp qtmp1;
          
          VectorXp &q = a.coeff[scaleIndex];

          if(gwt->getPsi().N() != 0){
            qtmp1 = gwt->getPsi() * q.head(gwt->getPsi().cols()); 
            qtmp1 += gwt->getW();
            
          }
          else{
            qtmp1 = gwt->getW();
          }
          
          VectorXp qtmp2 = parent->getPhi() * ( parent->getPhi().transpose() * sumq);
       
          qtmp1 -= qtmp2;
          sumq += qtmp1; 
          
          gwt = parent;
          parent = dynamic_cast< GWTNode<TPrecision> *>( gwt->getParent() );
          scaleIndex--;
        }
      }
      else{
        sumq *=  0;
      }

      VectorXp &p = a.coeff[scaleIndex];
      VectorXp x = (gwt->getPhi() * p ) + gwt->getCenter() + sumq;
      

      return x;
    };



};

#endif
