#ifndef WASSERSTEINNODEDISTANCE_H
#define WASSERSTEINNODEDISTANCE_H


#include "NodeDistance.h"
#include "GWT.h"
#include "EigenWasserstein.h"





template < typename TPrecision >
class WassersteinNodeDistance : public NodeDistance<TPrecision> {
  
  public:  
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> MatrixXp;
    typedef typename Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> VectorXp;

  private:
    Wasserstein<TPrecision> wstein;

  public:

    WassersteinNodeDistance(){};
    ~WassersteinNodeDistance(){};

    TPrecision distance(GMRANode<TPrecision> *n1, GMRANode<TPrecision> *n2){
      GWTNode<TPrecision> *gwt1 = dynamic_cast<GWTNode<TPrecision> *>( n1 );
      if(gwt1 == NULL){
        GMRANodeDecorator<TPrecision> *dec =
          dynamic_cast< GMRANodeDecorator<TPrecision> *>( n1 );
        while(dec != NULL && gwt1==NULL){
           n1 = dec->getDecoratedNode();
           gwt1 = dynamic_cast<GWTNode<TPrecision> *>( n1 ); 
           dec = dynamic_cast< GMRANodeDecorator<TPrecision> *>( n1 );

        }
      }

      GWTNode<TPrecision> *gwt2 = dynamic_cast<GWTNode<TPrecision> *>( n2 );
      if(gwt2 == NULL){
        GMRANodeDecorator<TPrecision> *dec =
          dynamic_cast< GMRANodeDecorator<TPrecision> *>( n2 );
        while(dec != NULL && gwt2==NULL){
          n2 = dec->getDecoratedNode();
          gwt2 = dynamic_cast<GWTNode<TPrecision> *>( n2 ); 
          dec = dynamic_cast< GMRANodeDecorator<TPrecision> *>( n2 );

        }
      }


      VectorXp sigma12 = gwt1->getSigma().array().square();
      VectorXp sigma22 = gwt2->getSigma().array().square();


      TPrecision d = wstein.distance( gwt1->getPhi(), sigma12,
          gwt1->getCenter(), gwt2->getPhi(), sigma22, gwt2->getCenter() );
      return d;
    };

};




#endif

