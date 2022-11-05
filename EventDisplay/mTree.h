//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 27 10:30:31 2022 by ROOT version 6.22/09
// from TTree mTree/Track Tree
// found on file: test.root
//////////////////////////////////////////////////////////

#ifndef mTree_h
#define mTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class mTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           mRunId;
   Int_t           mEvtId;
   Int_t           mEvtIndex;
   Float_t         mVx;
   Float_t         mVy;
   Float_t         mVz;
   Short_t         mRefMult;
   Int_t           mNT;
   Short_t         mNHits[10000];   //[mNT]
   Float_t         mGPx[10000];   //[mNT]
   Float_t         mGPy[10000];   //[mNT]
   Float_t         mGPz[10000];   //[mNT]
   Float_t         mOx[10000];   //[mNT]
   Float_t         mOy[10000];   //[mNT]
   Float_t         mOz[10000];   //[mNT]
   Float_t         mBTofX[10000];   //[mNT]
   Float_t         mBTofY[10000];   //[mNT]
   Float_t         mBTofZ[10000];   //[mNT]

   // List of branches
   TBranch        *b_mRunId;   //!
   TBranch        *b_mEvtId;   //!
   TBranch        *b_mEvtIndex;   //!
   TBranch        *b_mVx;   //!
   TBranch        *b_mVy;   //!
   TBranch        *b_mVz;   //!
   TBranch        *b_mRefMult;   //!
   TBranch        *b_mNT;   //!
   TBranch        *b_mNHits;   //!
   TBranch        *b_mGPx;   //!
   TBranch        *b_mGPy;   //!
   TBranch        *b_mGPz;   //!
   TBranch        *b_mOx;   //!
   TBranch        *b_mOy;   //!
   TBranch        *b_mOz;   //!
   TBranch        *b_mBTofX;   //!
   TBranch        *b_mBTofY;   //!
   TBranch        *b_mBTofZ;   //!

   mTree(TTree *tree=0);
   virtual ~mTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

//#ifdef mTree_cxx
mTree::mTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("test.root");
      }
      f->GetObject("mTree",tree);

   }
   Init(tree);
}

mTree::~mTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t mTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t mTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void mTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("mRunId", &mRunId, &b_mRunId);
   fChain->SetBranchAddress("mEvtId", &mEvtId, &b_mEvtId);
   fChain->SetBranchAddress("mEvtIndex", &mEvtIndex, &b_mEvtIndex);
   fChain->SetBranchAddress("mVx", &mVx, &b_mVx);
   fChain->SetBranchAddress("mVy", &mVy, &b_mVy);
   fChain->SetBranchAddress("mVz", &mVz, &b_mVz);
   fChain->SetBranchAddress("mRefMult", &mRefMult, &b_mRefMult);
   fChain->SetBranchAddress("mNT", &mNT, &b_mNT);
   fChain->SetBranchAddress("mNHits", mNHits, &b_mNHits);
   fChain->SetBranchAddress("mGPx", mGPx, &b_mGPx);
   fChain->SetBranchAddress("mGPy", mGPy, &b_mGPy);
   fChain->SetBranchAddress("mGPz", mGPz, &b_mGPz);
   fChain->SetBranchAddress("mOx", mOx, &b_mOx);
   fChain->SetBranchAddress("mOy", mOy, &b_mOy);
   fChain->SetBranchAddress("mOz", mOz, &b_mOz);
   fChain->SetBranchAddress("mBTofX", mBTofX, &b_mBTofX);
   fChain->SetBranchAddress("mBTofY", mBTofY, &b_mBTofY);
   fChain->SetBranchAddress("mBTofZ", mBTofZ, &b_mBTofZ);
   Notify();
}

Bool_t mTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void mTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t mTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
//#endif // #ifdef mTree_cxx
