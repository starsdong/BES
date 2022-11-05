//
#include <stdio>
#include <iomanip>
#include "mTree.h"

class mTree;
mTree *t = 0;
//==============================================================================
// Constants.
//------------------------------------------------------------------------------
const Double_t kR_MAX = 250.;
const Double_t kZ_MAX = 250.;

const Double_t kR_TPC_min = 46;
const Double_t kR_TPC_max = 200;
const Double_t kZ_TPC_d   = 200;

const Double_t kR_BP_min = 2.00;
const Double_t kR_BP_max = 2.076;
const Double_t kZ_BP_d   = 30;

// Solenoid field along z, in Tesla.
const Double_t kMagField = 0.498948;

// Draw Options
const Color_t  kColors[3] = { kRed, kGreen, kYellow };
const Double_t kTrackHitOffset = 0.5;


//==============================================================================
// Global variables.
//------------------------------------------------------------------------------

// Implemented in MultiView.C
class MultiView;
MultiView* gMultiView = 0;

TEveTrackList *gTrackList = 0;
TEvePointSet  *gVtxList = 0;
TEvePointSet  *gHitList = 0;
TEvePointSet  *gTrackHitList = 0;
TRandom3 *gRandom = new TRandom3;
Int_t index = 0;
Int_t NMAX = 0;
TGeoManager *gGeoManager = new TGeoManager;

//==============================================================================
// Forward decalarations of CINT functions.
//------------------------------------------------------------------------------


//==============================================================================
// Main - evt_display()
//------------------------------------------------------------------------------

void evt_display(const double B = 1)
//, const Int_t runnumber = 15081042)
{
  if (gROOT->LoadMacro("$ROOTSYS/tutorials/eve/MultiView.C+") != 0)
    {
      Error("evt_display()", "Failed loading MultiView.C in compiled mode.");
      return;
    }

  //========================================================================
  //========================================================================
  TChain *chain = new TChain("mTree");
  char inname[100];
  //  sprintf(inname, "output/Event_%d.root",runnumber);
  //  chain->AddFile(inname);
  chain->AddFile("test.root");
  NMAX = chain->GetEntries();
  cout << " Total number of events = " << NMAX << endl;
  index = 0;
  
  t = new mTree(chain);
  
  //========================================================================
  // Create views and containers.
  //========================================================================
  
  TEveManager::Create();
  
  make_geometry();
  
  init(B);
  
  //========================================================================
  //========================================================================
  
  evt_make_gui();
  process_event(index);
  
  gEve->Redraw3D(kTRUE);
  
}

//==============================================================================
// intitalize track/hit lists
//------------------------------------------------------------------------------


void init(const double B)
{
  gEve->GetBrowser()->GetTabRight()->SetTab(1);
  gTrackList = new TEveTrackList("Rec Tracks"); 
  gTrackList->SetMainColor(kYellow);
  gTrackList->SetMarkerColor(kRed);
  gTrackList->SetMarkerStyle(4);
  gTrackList->SetMarkerSize(0.5);
  gEve->AddElement(gTrackList);
  
  TEveTrackPropagator* trkProp = gTrackList->GetPropagator();
  trkProp->SetMagField(kMagField*B);
  // trkProp->SetMaxR(kR_TPC_max);
  // trkProp->SetMaxZ(kZ_TPC_d);
  trkProp->SetMaxR(kR_MAX);
  trkProp->SetMaxZ(kZ_MAX);
  
  gVtxList = new TEvePointSet("Primary Vertex");
  gVtxList->SetMainColor(kRed);
  gVtxList->SetMarkerColor(kYellow);
  gVtxList->SetMarkerStyle(20);
  gVtxList->SetMarkerSize(1.0);
  
  // gHitList = new TEvePointSet("HFT Rec Hits");
  // gHitList->SetMainColor(kRed);
  // gHitList->SetMarkerColor(kWhite);
  // gHitList->SetMarkerStyle(20);
  // gHitList->SetMarkerSize(1.0);
  
  gTrackHitList = new TEvePointSet("Track Rec Hits");
  gTrackHitList->SetMainColor(kRed);
  gTrackHitList->SetMarkerColor(kYellow);
  gTrackHitList->SetMarkerStyle(20);
  gTrackHitList->SetMarkerSize(0.8);
  
}

//==============================================================================
// Next event
//------------------------------------------------------------------------------

void process_event(Int_t iEvt)
{
  if(iEvt>=NMAX) {
    cout << " End of the tree! Go backward! " << endl;
  } else if(iEvt<0) {
    cout << " Beginning of the tree! Go forward! " << endl;
  }
  
  cout << "begin " << index << "th entry...." << endl;
  t->GetEntry(iEvt);
  
  gTrackList->DestroyElements();
  gVtxList->Reset();
  //  gHitList->Reset();
  gTrackHitList->Reset();
  
  int runId = t->mRunId;
  int evtId = t->mEvtId;
  int evtIndex = t->mEvtIndex;
  
  // Load verteice/hits
  float vx = t->mVx;
  float vy = t->mVy;
  float vz = t->mVz;
  cout << "==> Event Index " << evtIndex << "\t RunId/EvtId = " << runId << "/" << evtId << "\t Vtx = " << vx << " " << vy << " " << vz << endl;
  //  gVtxList->SetNextPoint(t->fVertex[0], t->fVertex[1], t->fVertex[2]);
  //  gEve->AddElement(gVtxList);

  /*  
  Int_t nHits = t->fHits_;
  for(int j=0;j<nHits;j++) {
    int id = t->fHits_Id[j];
    int nRawHits = t->fHits_nRawHits[j];
    if(id<1000 && ( nRawHits <= kPXL_Cluster_Min ||  nRawHits > kPXL_Cluster_Max ) ) continue;
    
    if(id<1000 && Status_PXL[id-1]) continue; // remove noisy channels
    if(id>1000 && Status_IST[id-1-1000]) continue;

    //    cout << " Adding a new hit " << id << " " <<  t->fHits_xG[j] << " " << t->fHits_yG[j] << " " << t->fHits_zG[j] << endl;
    gHitList->SetNextPoint(t->fHits_xG[j], t->fHits_yG[j], t->fHits_zG[j]);
  }
  gEve->AddElement(gHitList);
  */
  
  
  // Load tracks
  TEveTrackPropagator *trkProp = gTrackList->GetPropagator();
  Int_t nTracks = t->mNT;
  cout << "    # of Tracks = " << nTracks << endl;
  for (Int_t j = 0; j < nTracks; ++j) {
    float ox = t->mOx[j];
    float oy = t->mOy[j];
    float oz = t->mOz[j];
    if(sqrt(ox*ox+oy*oy)>5.0) continue;
    if(fabs(oz)>100.0) continue;
    TEveVectorT<double> origin(t->mOx[j], t->mOy[j], t->mOz[j]);
    TEveVectorT<double> mom(t->mGPx[j], t->mGPy[j], t->mGPz[j]);
    Int_t charge = (t->mNHits[j]>0) ? +1 : -1;
    
    TEveRecTrackT<double> tR;
    tR.fIndex = j;
    tR.fP = mom;
    tR.fV = origin;
    tR.fSign = charge;
    
    TEveTrack* track = new TEveTrack(&tR, trkProp);
    track->SetName(Form("%s [%d]", "rec", j));
    track->SetStdTitle();
    track->SetAttLineAttMarker(gTrackList);
    if (charge == +1)
      track->SetLineColor(kColors[0]);
    else
      track->SetLineColor(kColors[1]);
    
    //    cout << " Adding a new track " << t->fTracks_fPx[j] << " " << t->fTracks_fPy[j] << " " << t->fTracks_fPz[j] << endl;
    gTrackList->AddElement(track);
        
    int nhits = abs(t->mNHits[j]);
    int nhits_tpc = nhits%100;    
    //    cout << " Number of hits on PXL1/PXL2/IST/SSD/TPC = " << nhits_pxl1 << "/" << nhits_pxl2 << "/" << nhits_ist << "/" << nhits_ssd << "/" << nhits_tpc << endl;
  }  
  gTrackList->MakeTracks();
  gEve->AddElement(gTrackHitList);
  
  gEve->SetStatusLine(Form("run#%d event#%d",runId,evtId));
  
  TEveElement* top = gEve->GetCurrentEvent();
  
  gMultiView->DestroyEventRPhi();
  gMultiView->ImportEventRPhi(top);
  
  gMultiView->DestroyEventRhoZ();
  gMultiView->ImportEventRhoZ(top);
  
  gEve->Redraw3D();
  
}

void selectDaughterVisible(TGeoNode *node, const char *name)
{
  int nn=node->GetVolume()->GetNdaughters();
  for(int i = 0; i< nn; i++) {
    TGeoNode *daughter = node->GetVolume()->GetNode(i);
    if(!daughter) continue;
    if(strstr(daughter->GetName(), name)!=0) {
      cout << "  Found this node " << daughter->GetName() << " set to be visible" << endl;
      daughter->GetVolume()->SetVisibility(1);
    } else {
      daughter->GetVolume()->SetVisibility(0);
      if(daughter->GetVolume()->GetNdaughters()!=0) {
	selectDaughterVisible(daughter, name);
      }
    }
  }
}

void make_geometry()
{
  TEveElementList *STAR = new TEveElementList("Geometry");
  
  // gROOT->LoadMacro("y2014.C");
  // y2014();
  gGeoManager = gEve->GetGeometry("files/star_2014.root");
  
  TGeoVolume* top = gGeoManager->GetTopVolume()->FindNode("CAVE_1")->GetVolume();
  
  TGeoNode *tpc_mom = top->FindNode("TPCE_1");
//  selectDaughterVisible(tpc_mom, "TPCM");  // Central Membrane
//  selectDaughterVisible(tpc_mom, "TSAW");  // 
//  selectDaughterVisible(tpc_mom, "TWMR");
//  selectDaughterVisible(tpc_mom, "TWRB"); // Sector ribs
  // selectDaughterVisible(tpc_mom, "TINX");  // Central Membrane
  // selectDaughterVisible(tpc_mom, "TONX");  // Central Membrane
  TEveGeoTopNode* tpc = new TEveGeoTopNode(gGeoManager, tpc_mom);
//  tpc->SetMainTransparency(80);
  tpc->SetVisLevel(5);
  STAR->AddElement(tpc);

  
  gMultiView = new MultiView;
  gMultiView->ImportGeomRPhi(STAR);
  gMultiView->ImportGeomRhoZ(STAR);
}
//==============================================================================
// GUI stuff
//------------------------------------------------------------------------------
class EvNavHandler
{
public:
  void Fwd()
  {
    process_event(++index);
  }
  void Bck()
  {
    index--;
    process_event(index);
  }
};

//______________________________________________________________________________
void evt_make_gui()
{
  // Create minimal GUI for event navigation.
  
  TEveBrowser* browser = gEve->GetBrowser();
  browser->StartEmbedding(TRootBrowser::kLeft);
  
  TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 400, 400);
  frmMain->SetWindowName("XX GUI");
  frmMain->SetCleanup(kDeepCleanup);
  
  TGHorizontalFrame* hf = new TGHorizontalFrame(frmMain);
  {
    TString icondir( Form("%s/icons/", gSystem->Getenv("ROOTSYS")) );
    TGPictureButton* b = 0;
    EvNavHandler    *fh = new EvNavHandler;
    
    b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoBack.gif"));
    //      b->SetEnabled(kFALSE);
    b->SetToolTipText("Go to previous event - not supported.");
    hf->AddFrame(b);
    b->Connect("Clicked()", "EvNavHandler", fh, "Bck()");
    
    b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoForward.gif"));
    b->SetToolTipText("Generate new event.");
    hf->AddFrame(b);
    b->Connect("Clicked()", "EvNavHandler", fh, "Fwd()");

  }
  frmMain->AddFrame(hf);
  
  frmMain->MapSubwindows();
  frmMain->Resize();
  frmMain->MapWindow();
  
  browser->StopEmbedding();
  browser->SetTabTitle("Event Control", 0);
}
