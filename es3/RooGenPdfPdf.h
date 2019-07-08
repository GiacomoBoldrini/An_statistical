/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOGENPDFPDF
#define ROOGENPDFPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooGenPdfPdf : public RooAbsPdf {
public:
  RooGenPdfPdf() {} ; 
  RooGenPdfPdf(const char *name, const char *title,
	      RooAbsReal& _t,
	      RooAbsReal& _p,
	      RooAbsReal& _alpha,
	      RooAbsReal& _beta,
	      RooAbsReal& _gamma);
  RooGenPdfPdf(const RooGenPdfPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooGenPdfPdf(*this,newname); }
  inline virtual ~RooGenPdfPdf() { }

protected:

  RooRealProxy t ;
  RooRealProxy p ;
  RooRealProxy alpha ;
  RooRealProxy beta ;
  RooRealProxy gamma ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooGenPdfPdf,1) // Your description goes here...
};
 
#endif