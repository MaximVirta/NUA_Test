// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME nanoDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "src/JBaseTrack.h"
#include "src/JBaseEventHeader.h"
#include "src/JBaseEventHeader.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_JBaseTrack(void *p = nullptr);
   static void *newArray_JBaseTrack(Long_t size, void *p);
   static void delete_JBaseTrack(void *p);
   static void deleteArray_JBaseTrack(void *p);
   static void destruct_JBaseTrack(void *p);
   static void streamer_JBaseTrack(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::JBaseTrack*)
   {
      ::JBaseTrack *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::JBaseTrack >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("JBaseTrack", ::JBaseTrack::Class_Version(), "src/JBaseTrack.h", 26,
                  typeid(::JBaseTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::JBaseTrack::Dictionary, isa_proxy, 16,
                  sizeof(::JBaseTrack) );
      instance.SetNew(&new_JBaseTrack);
      instance.SetNewArray(&newArray_JBaseTrack);
      instance.SetDelete(&delete_JBaseTrack);
      instance.SetDeleteArray(&deleteArray_JBaseTrack);
      instance.SetDestructor(&destruct_JBaseTrack);
      instance.SetStreamerFunc(&streamer_JBaseTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::JBaseTrack*)
   {
      return GenerateInitInstanceLocal((::JBaseTrack*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::JBaseTrack*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_JBaseEventHeader(void *p = nullptr);
   static void *newArray_JBaseEventHeader(Long_t size, void *p);
   static void delete_JBaseEventHeader(void *p);
   static void deleteArray_JBaseEventHeader(void *p);
   static void destruct_JBaseEventHeader(void *p);
   static void streamer_JBaseEventHeader(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::JBaseEventHeader*)
   {
      ::JBaseEventHeader *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::JBaseEventHeader >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("JBaseEventHeader", ::JBaseEventHeader::Class_Version(), "src/JBaseEventHeader.h", 22,
                  typeid(::JBaseEventHeader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::JBaseEventHeader::Dictionary, isa_proxy, 16,
                  sizeof(::JBaseEventHeader) );
      instance.SetNew(&new_JBaseEventHeader);
      instance.SetNewArray(&newArray_JBaseEventHeader);
      instance.SetDelete(&delete_JBaseEventHeader);
      instance.SetDeleteArray(&deleteArray_JBaseEventHeader);
      instance.SetDestructor(&destruct_JBaseEventHeader);
      instance.SetStreamerFunc(&streamer_JBaseEventHeader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::JBaseEventHeader*)
   {
      return GenerateInitInstanceLocal((::JBaseEventHeader*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::JBaseEventHeader*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr JBaseTrack::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *JBaseTrack::Class_Name()
{
   return "JBaseTrack";
}

//______________________________________________________________________________
const char *JBaseTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::JBaseTrack*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int JBaseTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::JBaseTrack*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *JBaseTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::JBaseTrack*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *JBaseTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::JBaseTrack*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr JBaseEventHeader::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *JBaseEventHeader::Class_Name()
{
   return "JBaseEventHeader";
}

//______________________________________________________________________________
const char *JBaseEventHeader::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::JBaseEventHeader*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int JBaseEventHeader::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::JBaseEventHeader*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *JBaseEventHeader::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::JBaseEventHeader*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *JBaseEventHeader::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::JBaseEventHeader*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void JBaseTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class JBaseTrack.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TLorentzVector::Streamer(R__b);
      R__b >> fID;
      R__b >> fLabel;
      R__b >> fParticleType;
      R__b >> fCharge;
      R__b >> fStatus;
      R__b.CheckByteCount(R__s, R__c, JBaseTrack::IsA());
   } else {
      R__c = R__b.WriteVersion(JBaseTrack::IsA(), kTRUE);
      TLorentzVector::Streamer(R__b);
      R__b << fID;
      R__b << fLabel;
      R__b << fParticleType;
      R__b << fCharge;
      R__b << fStatus;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_JBaseTrack(void *p) {
      return  p ? new(p) ::JBaseTrack : new ::JBaseTrack;
   }
   static void *newArray_JBaseTrack(Long_t nElements, void *p) {
      return p ? new(p) ::JBaseTrack[nElements] : new ::JBaseTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_JBaseTrack(void *p) {
      delete ((::JBaseTrack*)p);
   }
   static void deleteArray_JBaseTrack(void *p) {
      delete [] ((::JBaseTrack*)p);
   }
   static void destruct_JBaseTrack(void *p) {
      typedef ::JBaseTrack current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_JBaseTrack(TBuffer &buf, void *obj) {
      ((::JBaseTrack*)obj)->::JBaseTrack::Streamer(buf);
   }
} // end of namespace ROOT for class ::JBaseTrack

//______________________________________________________________________________
void JBaseEventHeader::Streamer(TBuffer &R__b)
{
   // Stream an object of class JBaseEventHeader.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      R__b >> fEventID;
      {float R_Dummy; R__b >> R_Dummy; fCentrality=Double32_t(R_Dummy);}
      {float R_Dummy; R__b >> R_Dummy; fEP=Double32_t(R_Dummy);}
      {float R_Dummy; R__b >> R_Dummy; fVtxX=Double32_t(R_Dummy);}
      {float R_Dummy; R__b >> R_Dummy; fVtxY=Double32_t(R_Dummy);}
      {float R_Dummy; R__b >> R_Dummy; fVtxZ=Double32_t(R_Dummy);}
      {float R_Dummy; R__b >> R_Dummy; fVtxZErr=Double32_t(R_Dummy);}
      R__b.CheckByteCount(R__s, R__c, JBaseEventHeader::IsA());
   } else {
      R__c = R__b.WriteVersion(JBaseEventHeader::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      R__b << fEventID;
      R__b << float(fCentrality);
      R__b << float(fEP);
      R__b << float(fVtxX);
      R__b << float(fVtxY);
      R__b << float(fVtxZ);
      R__b << float(fVtxZErr);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_JBaseEventHeader(void *p) {
      return  p ? new(p) ::JBaseEventHeader : new ::JBaseEventHeader;
   }
   static void *newArray_JBaseEventHeader(Long_t nElements, void *p) {
      return p ? new(p) ::JBaseEventHeader[nElements] : new ::JBaseEventHeader[nElements];
   }
   // Wrapper around operator delete
   static void delete_JBaseEventHeader(void *p) {
      delete ((::JBaseEventHeader*)p);
   }
   static void deleteArray_JBaseEventHeader(void *p) {
      delete [] ((::JBaseEventHeader*)p);
   }
   static void destruct_JBaseEventHeader(void *p) {
      typedef ::JBaseEventHeader current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_JBaseEventHeader(TBuffer &buf, void *obj) {
      ((::JBaseEventHeader*)obj)->::JBaseEventHeader::Streamer(buf);
   }
} // end of namespace ROOT for class ::JBaseEventHeader

namespace {
  void TriggerDictionaryInitialization_nanoDict_Impl() {
    static const char* headers[] = {
"src/JBaseTrack.h",
"src/JBaseEventHeader.h",
"src/JBaseEventHeader.h",
nullptr
    };
    static const char* includePaths[] = {
"/home/maxim/lib/root/include/",
"/home/maxim/Documents/Work/Machine_Learning/NUA_Test/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "nanoDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$src/JBaseTrack.h")))  JBaseTrack;
class __attribute__((annotate("$clingAutoload$src/JBaseEventHeader.h")))  JBaseEventHeader;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "nanoDict dictionary payload"

#ifndef one
  #define one 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "src/JBaseTrack.h"
#include "src/JBaseEventHeader.h"
#include "src/JBaseEventHeader.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"JBaseEventHeader", payloadCode, "@",
"JBaseTrack", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("nanoDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_nanoDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_nanoDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_nanoDict() {
  TriggerDictionaryInitialization_nanoDict_Impl();
}
