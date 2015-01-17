//
// File generated by rootcint at Fri Jan 16 18:16:48 2015

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME MCDimuonRecoDict
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "MCDimuonRecoDict.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"

// Direct notice to TROOT of the dictionary's loading.
namespace {
   static struct DictInit {
      DictInit() {
         ROOT::RegisterModule();
      }
   } __TheDictionaryInitializer;
}

// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void MCDimuonReco_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void delete_MCDimuonReco(void *p);
   static void deleteArray_MCDimuonReco(void *p);
   static void destruct_MCDimuonReco(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MCDimuonReco*)
   {
      ::MCDimuonReco *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MCDimuonReco >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MCDimuonReco", ::MCDimuonReco::Class_Version(), "./MCDimuonReco.h", 4,
                  typeid(::MCDimuonReco), DefineBehavior(ptr, ptr),
                  &::MCDimuonReco::Dictionary, isa_proxy, 4,
                  sizeof(::MCDimuonReco) );
      instance.SetDelete(&delete_MCDimuonReco);
      instance.SetDeleteArray(&deleteArray_MCDimuonReco);
      instance.SetDestructor(&destruct_MCDimuonReco);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MCDimuonReco*)
   {
      return GenerateInitInstanceLocal((::MCDimuonReco*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::MCDimuonReco*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *MCDimuonReco::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *MCDimuonReco::Class_Name()
{
   return "MCDimuonReco";
}

//______________________________________________________________________________
const char *MCDimuonReco::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MCDimuonReco*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MCDimuonReco::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MCDimuonReco*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void MCDimuonReco::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MCDimuonReco*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *MCDimuonReco::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MCDimuonReco*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void MCDimuonReco::Streamer(TBuffer &R__b)
{
   // Stream an object of class MCDimuonReco.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MCDimuonReco::Class(),this);
   } else {
      R__b.WriteClassBuffer(MCDimuonReco::Class(),this);
   }
}

//______________________________________________________________________________
void MCDimuonReco::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class MCDimuonReco.
      TClass *R__cl = ::MCDimuonReco::IsA();
      if (R__cl || R__insp.IsA()) { }
      R__insp.Inspect(R__cl, R__insp.GetParent(), "fmass", &fmass);
      TObject::ShowMembers(R__insp);
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_MCDimuonReco(void *p) {
      delete ((::MCDimuonReco*)p);
   }
   static void deleteArray_MCDimuonReco(void *p) {
      delete [] ((::MCDimuonReco*)p);
   }
   static void destruct_MCDimuonReco(void *p) {
      typedef ::MCDimuonReco current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MCDimuonReco

/********************************************************
* MCDimuonRecoDict.cxx
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableMCDimuonRecoDict();

extern "C" void G__set_cpp_environmentMCDimuonRecoDict() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("MCDimuonReco.h");
  G__cpp_reset_tagtableMCDimuonRecoDict();
}
#include <new>
extern "C" int G__cpp_dllrevMCDimuonRecoDict() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* MCDimuonReco */
static int G__MCDimuonRecoDict_169_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   MCDimuonReco* p = NULL;
   char* gvp = (char*) G__getgvp();
   //m: 1
   if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
     p = new MCDimuonReco(*(TreeHandler*) libp->para[0].ref);
   } else {
     p = new((void*) gvp) MCDimuonReco(*(TreeHandler*) libp->para[0].ref);
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__MCDimuonRecoDictLN_MCDimuonReco));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MCDimuonRecoDict_169_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 103, (long) ((MCDimuonReco*) G__getstructoffset())->findHighPtDimuon());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MCDimuonRecoDict_169_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 102, (double) ((const MCDimuonReco*) G__getstructoffset())->getDimuonMass());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MCDimuonRecoDict_169_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) MCDimuonReco::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MCDimuonRecoDict_169_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) MCDimuonReco::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MCDimuonRecoDict_169_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) MCDimuonReco::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MCDimuonRecoDict_169_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      MCDimuonReco::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MCDimuonRecoDict_169_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((MCDimuonReco*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MCDimuonRecoDict_169_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) MCDimuonReco::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MCDimuonRecoDict_169_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MCDimuonReco::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MCDimuonRecoDict_169_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) MCDimuonReco::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__MCDimuonRecoDict_169_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) MCDimuonReco::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__MCDimuonRecoDict_169_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   MCDimuonReco* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new MCDimuonReco(*(MCDimuonReco*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__MCDimuonRecoDictLN_MCDimuonReco));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef MCDimuonReco G__TMCDimuonReco;
static int G__MCDimuonRecoDict_169_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 1
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (MCDimuonReco*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((MCDimuonReco*) (soff+(sizeof(MCDimuonReco)*i)))->~G__TMCDimuonReco();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (MCDimuonReco*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((MCDimuonReco*) (soff))->~G__TMCDimuonReco();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* MCDimuonReco */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncMCDimuonRecoDict {
 public:
  G__Sizep2memfuncMCDimuonRecoDict(): p(&G__Sizep2memfuncMCDimuonRecoDict::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncMCDimuonRecoDict::*p)();
};

size_t G__get_sizep2memfuncMCDimuonRecoDict()
{
  G__Sizep2memfuncMCDimuonRecoDict a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceMCDimuonRecoDict() {

   /* Setting up class inheritance */
   if(0==G__getnumbaseclass(G__get_linked_tagnum(&G__MCDimuonRecoDictLN_MCDimuonReco))) {
     MCDimuonReco *G__Lderived;
     G__Lderived=(MCDimuonReco*)0x1000;
     {
       TObject *G__Lpbase=(TObject*)G__Lderived;
       G__inheritance_setup(G__get_linked_tagnum(&G__MCDimuonRecoDictLN_MCDimuonReco),G__get_linked_tagnum(&G__MCDimuonRecoDictLN_TObject),(long)G__Lpbase-(long)G__Lderived,1,1);
     }
   }
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableMCDimuonRecoDict() {

   /* Setting up typedef entry */
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__MCDimuonRecoDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__MCDimuonRecoDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__MCDimuonRecoDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__MCDimuonRecoDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__MCDimuonRecoDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__MCDimuonRecoDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__MCDimuonRecoDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__MCDimuonRecoDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__MCDimuonRecoDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__MCDimuonRecoDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* MCDimuonReco */
static void G__setup_memvarMCDimuonReco(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__MCDimuonRecoDictLN_MCDimuonReco));
   { MCDimuonReco *p; p=(MCDimuonReco*)0x1000; if (p) { }
   G__memvar_setup((void*)0,117,1,0,G__get_linked_tagnum(&G__MCDimuonRecoDictLN_TreeHandler),-1,-1,4,"fT=",0,(char*)NULL);
   G__memvar_setup((void*)0,102,0,0,-1,-1,-1,4,"fmass=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__MCDimuonRecoDictLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarMCDimuonRecoDict() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncMCDimuonReco(void) {
   /* MCDimuonReco */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__MCDimuonRecoDictLN_MCDimuonReco));
   G__memfunc_setup("MCDimuonReco",1157,G__MCDimuonRecoDict_169_0_1, 105, G__get_linked_tagnum(&G__MCDimuonRecoDictLN_MCDimuonReco), -1, 0, 1, 1, 1, 0, "u 'TreeHandler' - 1 - in", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("findHighPtDimuon",1617,G__MCDimuonRecoDict_169_0_2, 103, -1, -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("getDimuonMass",1344,G__MCDimuonRecoDict_169_0_3, 102, -1, -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__MCDimuonRecoDict_169_0_4, 85, G__get_linked_tagnum(&G__MCDimuonRecoDictLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&MCDimuonReco::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__MCDimuonRecoDict_169_0_5, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&MCDimuonReco::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__MCDimuonRecoDict_169_0_6, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&MCDimuonReco::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__MCDimuonRecoDict_169_0_7, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&MCDimuonReco::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,(G__InterfaceMethod) NULL,85, G__get_linked_tagnum(&G__MCDimuonRecoDictLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,(G__InterfaceMethod) NULL,121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__MCDimuonRecoDict_169_0_11, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__MCDimuonRecoDict_169_0_12, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&MCDimuonReco::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__MCDimuonRecoDict_169_0_13, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&MCDimuonReco::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__MCDimuonRecoDict_169_0_14, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&MCDimuonReco::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__MCDimuonRecoDict_169_0_15, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&MCDimuonReco::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("MCDimuonReco", 1157, G__MCDimuonRecoDict_169_0_16, (int) ('i'), G__get_linked_tagnum(&G__MCDimuonRecoDictLN_MCDimuonReco), -1, 0, 1, 1, 1, 0, "u 'MCDimuonReco' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~MCDimuonReco", 1283, G__MCDimuonRecoDict_169_0_17, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncMCDimuonRecoDict() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {
}

static void G__cpp_setup_global2() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalMCDimuonRecoDict() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
  G__cpp_setup_global2();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcMCDimuonRecoDict() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__MCDimuonRecoDictLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__MCDimuonRecoDictLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__MCDimuonRecoDictLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__MCDimuonRecoDictLN_TObject = { "TObject" , 99 , -1 };
G__linked_taginfo G__MCDimuonRecoDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__MCDimuonRecoDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__MCDimuonRecoDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__MCDimuonRecoDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__MCDimuonRecoDictLN_TreeHandler = { "TreeHandler" , 99 , -1 };
G__linked_taginfo G__MCDimuonRecoDictLN_MCDimuonReco = { "MCDimuonReco" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableMCDimuonRecoDict() {
  G__MCDimuonRecoDictLN_TClass.tagnum = -1 ;
  G__MCDimuonRecoDictLN_TBuffer.tagnum = -1 ;
  G__MCDimuonRecoDictLN_TMemberInspector.tagnum = -1 ;
  G__MCDimuonRecoDictLN_TObject.tagnum = -1 ;
  G__MCDimuonRecoDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__MCDimuonRecoDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__MCDimuonRecoDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__MCDimuonRecoDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__MCDimuonRecoDictLN_TreeHandler.tagnum = -1 ;
  G__MCDimuonRecoDictLN_MCDimuonReco.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableMCDimuonRecoDict() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__MCDimuonRecoDictLN_TClass);
   G__get_linked_tagnum_fwd(&G__MCDimuonRecoDictLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__MCDimuonRecoDictLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__MCDimuonRecoDictLN_TObject);
   G__get_linked_tagnum_fwd(&G__MCDimuonRecoDictLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__MCDimuonRecoDictLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__MCDimuonRecoDictLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__MCDimuonRecoDictLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__MCDimuonRecoDictLN_TreeHandler);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__MCDimuonRecoDictLN_MCDimuonReco),sizeof(MCDimuonReco),-1,323584,(char*)NULL,G__setup_memvarMCDimuonReco,G__setup_memfuncMCDimuonReco);
}
extern "C" void G__cpp_setupMCDimuonRecoDict(void) {
  G__check_setup_version(30051515,"G__cpp_setupMCDimuonRecoDict()");
  G__set_cpp_environmentMCDimuonRecoDict();
  G__cpp_setup_tagtableMCDimuonRecoDict();

  G__cpp_setup_inheritanceMCDimuonRecoDict();

  G__cpp_setup_typetableMCDimuonRecoDict();

  G__cpp_setup_memvarMCDimuonRecoDict();

  G__cpp_setup_memfuncMCDimuonRecoDict();
  G__cpp_setup_globalMCDimuonRecoDict();
  G__cpp_setup_funcMCDimuonRecoDict();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncMCDimuonRecoDict();
  return;
}
class G__cpp_setup_initMCDimuonRecoDict {
  public:
    G__cpp_setup_initMCDimuonRecoDict() { G__add_setup_func("MCDimuonRecoDict",(G__incsetup)(&G__cpp_setupMCDimuonRecoDict)); G__call_setup_funcs(); }
   ~G__cpp_setup_initMCDimuonRecoDict() { G__remove_setup_func("MCDimuonRecoDict"); }
};
G__cpp_setup_initMCDimuonRecoDict G__cpp_setup_initializerMCDimuonRecoDict;
