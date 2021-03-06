#ifndef STEPFLATTENER_H
#define STEPFLATTENER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Class to perform flattening of the event plane distribution
//       It stores necessary parameterizations and apply when requested
//
//*-- Author: Dmitri Peressounko (RRC KI)

// --- ROOT system ---
class TH2 ;
#include "TNamed.h"

class StEPFlattener : public TNamed {

 public:
  
  StEPFlattener() ;
  StEPFlattener(const char * name, Int_t v=2) ; //To separate different EP detectors use names
  StEPFlattener(const StEPFlattener & fl) ; 
  virtual ~StEPFlattener() ;
  StEPFlattener & operator = (const StEPFlattener & flat);

public:

  Double_t MakeFlat(Double_t oldPhi, Double_t centrality)const ; // Apply (centrality-dependent) flattening to oldPhi
  void SetParameterization(TH2 * h) ;  // Set Parameterization to use (see code for the meaning of parameters

private:
  Int_t fNCentrBins ; // Number of centrality bins
  Int_t fNHarmonics ; // Number of harmonics used in parameterization
  Int_t fNparam ;     // Total number of parameters (fNCentrBins*fNHarmonics)
  Int_t fV3 ;         // Use v2 or V3 flattening
  Double32_t *fParam ;  // [fNparam][-1.,1.,16] array of flattening parameters

  ClassDef(StEPFlattener,1) 

} ;

#endif //  STEPFLATTENER_H
