#ifndef GMRADECORATOR_H
#define GMRADECORATOR_H

#include "GMRANode.h"



//Decorator to add addiotnal functionality to an arbitrary GMRA tree
template <typename TPrecision>
class Decorator{
  public:
    //
    virtual GMRANode<TPrecision> *decorate(GMRANode<TPrecision> *node, GMRANode<TPrecision> *parent) = 0;

};









#endif
