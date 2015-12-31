#include "Precision.h"

#include <stdio.h>


#include "LinalgIO.h"
#include "IPCAGWT.h"
#include "DenseMatrix.h"
#include "CmdLine.h"

#include <set>

int main(int argc, char **argv){
  using namespace FortranLinalg;

  //Command line parsing
  TCLAP::CmdLine cmd("IPCA GWT", ' ', "1");

  TCLAP::ValueArg<std::string> tArg("t","tree", "Tree data file", true, "",
      "tree data file");
  cmd.add(tArg);

  TCLAP::ValueArg<std::string> outArg("o","out", "data file name", true, "", "file prefix to store results in"); 
  cmd.add(outArg);

  TCLAP::ValueArg<std::string> xArg("x","data", "matrix header file", false, "", "input data matrix"); 
  cmd.add(xArg);

  TCLAP::MultiArg<int> sArg("s","scales", "reconstruct at gmra scales (0 coarsest)", false, "integer"); 
  cmd.add(sArg);

  TCLAP::MultiArg<double> zArg("z","zero", "threshold for GWT coefficents", false, "double"); 
  cmd.add(zArg);

  TCLAP::SwitchArg dArg("d", "dict", "save full dictionary");
  cmd.add(dArg);

  TCLAP::SwitchArg lArg("l", "leaf", "save active leafs dictionary");
  cmd.add(lArg);

  try{ 
    cmd.parse( argc, argv );
  } 
  catch (TCLAP::ArgException &e){ 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    return -1;
  }


  std::string outFileName = outArg.getValue();
  std::string treeFileName = tArg.getValue();
  std::string xFile = xArg.getValue();
  std::vector<int> scales = sArg.getValue();
  std::vector<double> thresholds = zArg.getValue();
  std::sort( thresholds.begin(), thresholds.end() );


  IPCATree<Precision> ipcaTree;

  std::ifstream treeFile(treeFileName.c_str(), std::ios_base::in | std::ios_base::binary);
  ipcaTree.unflatten(treeFile);

  treeFile.close();

  std::cout << "Read tree " << std::endl;


  IPCAGWT<Precision> gwt;
  gwt.setTree(&ipcaTree);

  std::cout << "Constructed GWT information" << std::endl;


  if(scales.size() > 0 || thresholds.size() > 0){
    if(xFile.size() == 0){
      std::cerr << "No data supplied to do GWT transfrom of" << std::endl;
      exit(1);
    }


    DenseMatrix<Precision> X = LinalgIO<Precision>::readMatrix(xFile);
    DenseVector<Precision> x(X.M());
    DenseMatrix<Precision> Xr( X.M(), X.N()) ;
    std::vector< GWTCoefficients<Precision> > cs;
    int maxScale = 0;
    int d = 0;

    std::set<GMRANode<Precision> *> active;

    Precision maxCoeff = 0;
    for(int i=0; i<X.N(); i++){
      Linalg<Precision>::ExtractColumn(X, i, x);

      GWTCoefficients<Precision> c = gwt.encode(x);
      cs.push_back(c);

      for(int k=1; k< c.coeff.size(); k++){
        for(int l=0; l< c.coeff[k].N(); l++){
          if(maxCoeff < fabs(c.coeff[k](l))){
            maxCoeff = fabs(c.coeff[k](l));
          }
        }
      }
      if(maxScale < c.coeff.size()){
        maxScale = c.coeff.size();
      }
      if(d < c.maxD){
        d = c.maxD;
      } 
    }


    DenseMatrix<Precision> coeff(maxScale, X.N()*d);
    Linalg<Precision>::Zero(coeff);
    DenseMatrix<int> ids(maxScale, X.N());
    Linalg<int>::Set(ids, -1);

    for(int i=0; i<X.N(); i++){
      int index = i*d;
      GWTCoefficients<Precision> &c = cs[i];
      for(int m=0; m < c.coeff.size(); m++){
        ids(m, i) = c.ids[m];
        for(int n=0; n < c.coeff[m].N(); n++){
          coeff(m, index + n) = c.coeff[m](n);
        }
      }
    }


    std::stringstream ss2;
    ss2 << outFileName <<  "-gwt-coeff.data";
    LinalgIO<Precision>::writeMatrix(ss2.str(), coeff);

    std::stringstream ss3;
    ss3 << outFileName <<  "-gwt-coeff-ids.data";
    LinalgIO<int>::writeMatrix(ss3.str(), ids);

    // -------------- Scales -------------- //
    for(int k=0; k< scales.size(); k++){
      int scale = scales[k];

      std::set<GMRANode<Precision> *> active;

      for(int i=0; i<X.N(); i++){
        Linalg<Precision>::ExtractColumn(X, i, x);

        GWTCoefficients<Precision> &c = cs[i];

        DenseVector<Precision> xr = gwt.decode(c, scale);

        Linalg<Precision>::SetColumn(Xr, i, xr);
        xr.deallocate();
      }

      std::stringstream ss1;
      ss1 << outFileName << "-scale-" << scale << "-reconstructed.data";
      LinalgIO<Precision>::writeMatrix(ss1.str(), Xr);

    }




    // -------- thresholds -------- //
    for(int k=0; k < thresholds.size(); k++){
      double t = thresholds[k]*maxCoeff;
      std::cout << t << std::endl;

      for(int i=0; i<X.N(); i++){
        Linalg<Precision>::ExtractColumn(X, i, x);


        GWTCoefficients<Precision> &c = cs[i];

        for(int m=1; m<c.coeff.size(); m++){
          for(int n=0; n<c.coeff[m].N(); n++){
            if( fabs( c.coeff[m](n) ) < t ){
              c.coeff[m](n) = 0;
            }
          }
        }
        DenseVector<Precision> xr = gwt.decode(c);

        Linalg<Precision>::SetColumn(Xr, i, xr);
        xr.deallocate();

      }

      std::stringstream ss1;
      ss1 << outFileName << "-threshold-" << k << "-reconstructed.data";
      LinalgIO<Precision>::writeMatrix(ss1.str(), Xr);

    }


  }




  //Save all dictionary elements
  if( dArg.getValue() ){
    //store centers and phi for each scale
    class StoreDictionary : public Visitor<Precision>{
      private:
        int scale;
        typedef std::map<int, DenseVector<Precision> > CenterMap;
        CenterMap centers;
        typedef std::map<int, DenseMatrix<Precision> > PhiMap;
        PhiMap phis;
        typedef std::map<int, int > ScaleMap;
        ScaleMap scales;


      public:

        StoreDictionary() {
        };


        virtual void visit(GMRANode<Precision> *node){
          GWTNode<Precision> *info = (GWTNode<Precision> *) node;
          centers[info->getID()] = info->getCenter();
          phis[info->getID()] = info->getPhi();
          scales[info->getID()] = node->getScale();
        };

        DenseMatrix<Precision> getPhis(){
          int n = 0;
          for(PhiMap::iterator it=phis.begin(); it !=phis.end(); ++it){
            n += it->second.N();
          };

          DenseMatrix<Precision> M(phis[0].M(), n);
          int index = 0;
          for(PhiMap::iterator it=phis.begin(); it !=phis.end(); ++it){
           
            Linalg<Precision>::SetColumns(M, index, index+it->second.N(),
                it->second, 0);
            index += it->second.N();
          }
          return M;
        };

        DenseMatrix<Precision> getCenters(){
          DenseMatrix<Precision> M(centers[0].N(), centers.size());
          for(CenterMap::iterator it = centers.begin(); it != centers.end();
              ++it){
            Linalg<Precision>::SetColumn(M, it->first, it->second);
          }
          return M;
        };

        DenseVector<int> getDimensions(){
          DenseVector<int> dims(phis.size());
          for(PhiMap::iterator it=phis.begin(); it !=phis.end(); ++it){
            dims(it->first) = it->second.N();
          };
          return dims;
        };

        DenseVector<int> getScales(){
          DenseVector<int> scale(phis.size());
          for(ScaleMap::iterator it=scales.begin(); it !=scales.end(); ++it){
            scale(it->first) = it->second;
          };
          return scale;
        };

    };


    StoreDictionary storeThingyDo;
    ipcaTree.breadthFirstVisitor(&storeThingyDo);
    DenseMatrix<Precision> Dc = storeThingyDo.getCenters();
    DenseMatrix<Precision> Dphi = storeThingyDo.getPhis();
    DenseVector<int> Dd = storeThingyDo.getDimensions();
    DenseVector<int> Ds = storeThingyDo.getScales();
    //DenseVector<int> Dids = storeThingyDo.getIDs();

    std::stringstream ss1;
    ss1 << outFileName << "-dict-center.data";
    LinalgIO<Precision>::writeMatrix(ss1.str(), Dc );

    std::stringstream ss2;
    ss2 << outFileName << "-dict-phi.data";
    LinalgIO<Precision>::writeMatrix(ss2.str(), Dphi );

    std::stringstream ss3;
    ss3 << outFileName << "-dict-dims.data";
    LinalgIO<int>::writeVector(ss3.str(), Dd);

    std::stringstream ss4;
    ss4 << outFileName << "-dict-scales.data";
    LinalgIO<int>::writeVector(ss4.str(), Ds);

   
    //Ids are sorted from 0  to n in linear order 
    //std::stringstream ss5;
    //ss5 << outFileName << "-dict-idss.data";
    //LinalgIO<int>::writeVector(ss5.str(), Dids);

  }





  return 0;

}
