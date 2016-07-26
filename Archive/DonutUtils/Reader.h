// @(#)root/tmva $Id: Reader.h,v 1.34 2008/03/29 21:57:21 andreas.hoecker Exp $ 
// Author: Andreas Hoecker, Joerg Stelzer, Helge Voss, Kai Voss 

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : Reader                                                                *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Reader class to be used in the user application to interpret the trained  *
 *      MVAs in an analysis context                                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Andreas Hoecker <Andreas.Hocker@cern.ch> - CERN, Switzerland              *
 *      Xavier Prudent  <prudent@lapp.in2p3.fr>  - LAPP, France                   *
 *      Helge Voss      <Helge.Voss@cern.ch>     - MPI-K Heidelberg, Germany      *
 *      Kai Voss        <Kai.Voss@cern.ch>       - U. of Victoria, Canada         *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         * 
 *      U. of Victoria, Canada                                                    * 
 *      MPI-K Heidelberg, Germany                                                 * 
 *      LAPP, Annecy, France                                                      *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://ttmva.sourceforge.net/LICENSE)                                         *
 **********************************************************************************/

#ifndef ROOT_TMVA_Reader
#define ROOT_TMVA_Reader

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Reader                                                               //
//                                                                      //
// Reader class to be used in the user application to interpret the     //
// trained MVAs in an analysis context                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TMVA_Types
#include "TMVA/Types.h"
#endif
#ifndef ROOT_TMVA_IMethod
#include "TMVA/MethodBase.h"
#endif
#ifndef ROOT_TMVA_DataSet
#include "TMVA/DataSet.h"
#endif
#ifndef ROOT_TMVA_MsgLogger
#include "TMVA/MsgLogger.h"
#endif
#ifndef ROOT_TMVA_Configurable
#include "TMVA/Configurable.h"
#endif

#include <vector>
#include <map>
#include <string>

namespace TMVA {

   class Reader : public Configurable {

   public:
      
      // without prior specification of variables
      Reader( TString theOption="", Bool_t verbose = 0 );

      // STL types
      Reader( std::vector<std::string>&  varNames, TString theOption = "", Bool_t verbose = 0 );
      Reader( const std::string varNames, TString theOption, Bool_t verbose = 0 );  // format: "var1:var2:..."

      // Root types
      Reader( std::vector<TString>& varNames, TString theOption = "", Bool_t verbose = 0 );
      Reader( const TString varNames, TString theOption, Bool_t verbose = 0 );  // format: "var1:var2:..."

      virtual ~Reader( void );
  
      IMethod* BookMVA( const TString& methodTag, const TString& weightfile );
      IMethod* FindMVA( const TString& methodTag );

      Double_t EvaluateMVA( const std::vector<Float_t> &, const TString& methodTag, Double_t aux = 0 );    
      Double_t EvaluateMVA( const std::vector<Double_t>&, const TString& methodTag, Double_t aux = 0 );    
      Double_t EvaluateMVA( MethodBase* method,           Double_t aux = 0 );    
      Double_t EvaluateMVA( const TString& methodTag,            Double_t aux = 0 );    

      Double_t GetProba ( const TString& methodTag, Double_t ap_sig=0.5, Double_t mvaVal=-9999999 ); 
      Double_t GetRarity( const TString& methodTag, Double_t mvaVal=-9999999 );

      // accessors 
      virtual const char* GetName() const { return "Reader"; }
      Bool_t   Verbose( void ) const  { return fVerbose; }
      void     SetVerbose( Bool_t v ) { fVerbose = v; }

      const DataSet& Data() const { return *fDataSet; }
      DataSet& Data() { return *fDataSet; }
      
      // add and retrieve variables
      void AddVariable( const TString& expression, Float_t* );
      void AddVariable( const TString& expression, Int_t* );  
      const TString& GetVarName( Int_t ivar ) { return Data().GetExpression(ivar); }

   private:

      // this booking method is internal
      IMethod* BookMVA( Types::EMVA method,  TString weightfile );

      DataSet* fDataSet; // the data set
  
      // Init Reader class
      void Init( void );

      // Decode Constructor string (or TString) and fill variable name std::vector
      void DecodeVarNames( const std::string varNames );
      void DecodeVarNames( const TString varNames );

      // Reads method name and title from the weightfile
      void GetMethodNameTitle(const TString& weightfile, TString& methodName, TString& methodTitle);


      void DeclareOptions();

      Bool_t fVerbose;    // verbosity

      std::map<TString, IMethod*> fMethodMap; // map of methods

      mutable MsgLogger fLogger; // message logger

      ClassDef(Reader,0) // Interpret the trained MVAs in an analysis context
   };

}

#endif
