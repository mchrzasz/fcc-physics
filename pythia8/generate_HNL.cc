// Data model
#include "datamodel/EventInfo.h"
#include "datamodel/EventInfoCollection.h"
#include "datamodel/MCParticle.h"
#include "datamodel/MCParticleCollection.h"
#include "datamodel/GenJet.h"
#include "datamodel/GenJetCollection.h"
#include "datamodel/GenVertex.h"
#include "datamodel/GenVertexCollection.h"
#include "datamodel/LorentzVector.h"

// STL
#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>

// albers specific includes
#include "podio/EventStore.h"
#include "podio/ROOTWriter.h"

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include "Pythia8/ParticleData.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "TDatabasePDG.h"

// root stuff:
#include "TDatabasePDG.h"               // for TDatabasePDG
#include "TParticlePDG.h"


namespace CKM
{
      // a bit old fit, take the new one in the future
  static const double  Vud = 0.9742;
  static const double Vus = 0.2252;
  static const double    Vub = 4.2e-03;
  static const double    Vcd = 0.23;
  static const double    Vcs = 1.;
  static const double    Vcb = 4.1e-02;
  static const double    Vtd = 8.4e-03;
  static const double    Vts = 4.3e-02;
  static const double    Vtb = 0.89;


}
namespace constants
{
  static const double  GF=1.166379e-05;//GeV
  static const double  s2thetaw=0.23126;
  static const double pi=3.1415;

}

double decayConstant(int i)
{
  std::map<int,double> mymap;
  mymap[211]=0.1307;
  mymap[111]=0.130;
  mymap[113]=0.102;//rho0
  mymap[213]=0.102;  //rho+
  mymap[221]=1.2*0.130;  //eta
  mymap[331]=-0.45*0.130;  //eta(958)

  return mymap[i];

}
double CKMelemSq(int a)
{

  std::unordered_map<int, double> CKMelemSq = {{211,CKM::Vud*CKM::Vud}, //pi
                                               {113,CKM::Vud*CKM::Vud},  //rho
                                               {321,CKM::Vus*CKM::Vus}}  ; //K
  return CKMelemSq[a];
}

auto isEven = [](int x) -> bool
{
  return x % 2 == 0;
};






double CKMelemSq(int a, int b)
{
  //std::unordered_map< std::pair< int,int >, double> CKMelemSq =  { {{0,0},1.}};//, {{1,1},1.}, {{2,2},1.}, {{3,3},1.}, {{4,4},1.}, {{5,5},1. }  };
  // 0=u, 1=d, 2=s, 3=c, 4=b, 5=t
  if( a==b) return 1.;// the same flavour
  if( (a==0 && b==1) || (a==1 && b==0)) return CKM::Vud * CKM::Vud;
  if( (a==0 && b==2) || (a==2 && b==0)) return CKM::Vus * CKM::Vus;
  if( (a==0 && b==4) || (a==4 && b==0)) return CKM::Vub * CKM::Vub;

  if( (a==1 && b==3) || (a==3 && b==1)) return CKM::Vcd * CKM::Vcd;
  if( (a==2 && b==3) || (a==3 && b==2)) return CKM::Vcs * CKM::Vcs;
  if( (a==4 && b==3) || (a==3 && b==4)) return CKM::Vcb * CKM::Vcb;

  if( (a==1 && b==5) || (a==5 && b==1)) return CKM::Vtd * CKM::Vtd;
  if( (a==2 && b==5) || (a==5 && b==2)) return CKM::Vts * CKM::Vts;
  if( (a==4 && b==5) || (a==5 && b==4)) return CKM::Vtb * CKM::Vtb;

  else return 0.;

}















class couplings
{
public:
  double U2e, U2mu, U2tau;
  double u[3];
  couplings()
  {
    U2e=0.;
    U2mu=0.;
    U2tau=0.;
    u[0]= U2e;
    u[1]= U2mu;
    u[2]=U2tau;
  }
  couplings(double u2e, double u2mu, double u2tau)
  {
    U2e=u2e;
    U2mu=u2mu;
    U2tau=u2tau;
    u[0]= U2e;
    u[1]= U2mu;
    u[2]=U2tau;
  };
  double &operator[](int i)
  {
    return u[i];
  }


};







using namespace std;




class HNL
{
  // class thar provides calculations of Branching fractions and width
  // for different couplings and masses of HNL
  // M.Chrzaszcz, 09.04.2015
  HNL(couplings u, double Mass)
  {
    U2=u;
    U=couplings(sqrt(u.U2e), sqrt(u.U2mu), sqrt(u.U2tau));
    pdgBase = TDatabasePDG::Instance();
    m=Mass;


  };
  double Width_H_l(int h, int l)
  {
    //###################################################################
    // returns the decay width of HNL to charge meson and lepton
    // is the hadron int, int is lepton
    //###################################################################
    double alpha;
    if(abs(l)==11) alpha=1;
    if(abs(l)==13) alpha=2;
    if(abs(l)==15) alpha=3;

    TParticlePDG *tPartH = pdgBase->GetParticle(h);
    double h_m=tPartH->Mass();
    TParticlePDG *tPartl = pdgBase->GetParticle(l);
    double l_m=tPartl->Mass();

    if(m<h_m+l_m) return 0.; // no decay

    double width = (abs(U2[alpha-1])/(16.*constants::pi))*(constants::GF * constants::GF)*(decayConstant(h)*decayConstant(h));
    width = width*(pow(m,3.))*CKMelemSq(h);
    double par = ( ((1 - pow(l_m,2)) /(pow(m,2)  ))*((1 - pow(l_m,2)) /(pow(m,2)  ))
            - ( ( pow(h_m,2) )/( pow(m,2))* (1 + (( pow(l_m,2) )/(  pow(m,2) ))) ));
    if(h==113 || h==213)
    {
      par=( ((1 - ((pow(l_m,2))/(pow(m,2))))*(1 - ((pow(l_m,2))/(pow(m,2)))))
            + ( (pow(h_m,2))/(pow(m,2))
                * (1 + (((pow(l_m,2) - 2.*pow(h_m,2)))/(pow(m,2)))) ) );
      par = par*2./(pow(h_m,2));

    }
    width = width*par;
    double rad = sqrt( ( 1-((h_m-l_m)*(h_m-l_m))/(pow(m,2)) )
                       * ( ( 1-((h_m+l_m)*(h_m+l_m))/(pow(m,2)) ) ) );

    width = width*rad;
    width = 2.*width; // Majorana case (charge conjugate channels)
    return width;
    


  }
  double Width_3nu()
  {
    double width = (constants::GF* constants::GF)*(pow(m,5))*(U2[0]+U2[1]+U2[2])/(192.*( pow(constants::pi,3) ));
    width = 2.*width; //
    return width;
    
  }
  double Width_H0_nu(int h, int l)
  {
    TParticlePDG *tPartH = pdgBase->GetParticle(h);
    double h_m=tPartH->Mass();
    TParticlePDG *tPartl = pdgBase->GetParticle(l);
    double l_m=tPartl->Mass();
    
    double alpha;
    if(abs(l)==12) alpha=1;
    if(abs(l)==14) alpha=2;
    if(abs(l)==16) alpha=3;
    


    if( h_m < m) return 0.;
    double width = (abs(U2[alpha-1])/(32.*constants::pi))*(constants::GF*constants::GF )*(decayConstant(h)* decayConstant(h));
    double  par = (pow(m,3))*((  pow( (1 - pow(h_m,2)  /(pow(m,2))),2)));
    if(h==113 || h==213)
    {
      par = par*2./(pow(h_m,2));
      par = par*(1 + 2.*(pow(h_m,2))/(pow(m,2)));
            
    }
    
    
    
    width = width*par;
    width = 2.*width;
    return width;
        
  }
  

private:

  couplings U2;
  couplings U;
  double m;

  TDatabasePDG* pdgBase ;//= TDatabasePDG::Instance();




};








int main(int argc, char** argv) {
  if( argc != 2) {
    std::cerr<<"Usage: pythiafcc-generate <pythia card file>"<<std::endl;
    return 1;
  }
  const char* card_file = argv[1];
  std::string output(card_file);
  size_t dot = output.find_last_of(".");
  if(dot != std::string::npos){
    output.replace(dot, output.size()-dot, ".root");
  }
  else {
    output += ".root";
  }
  size_t slash = output.find_last_of("/");
  if(slash != std::string::npos){
    output.replace(0, slash+1, "");
  }
  std::cout<<"start processing2"<<std::endl;
  std::cout<<"output file: "<<output<<std::endl;

  auto store  = podio::EventStore();
  auto writer = podio::ROOTWriter(output, &store);


  fcc::EventInfoCollection& evinfocoll = store.create<fcc::EventInfoCollection>("EventInfo");
  fcc::MCParticleCollection& pcoll = store.create<fcc::MCParticleCollection>("GenParticle");
  fcc::GenVertexCollection& vcoll = store.create<fcc::GenVertexCollection>("GenVertex");
  fcc::GenJetCollection& genjetcoll = store.create<fcc::GenJetCollection>("GenJet");

  writer.registerForWrite<fcc::EventInfoCollection>("EventInfo");
  writer.registerForWrite<fcc::MCParticleCollection>("GenParticle");
  writer.registerForWrite<fcc::GenVertexCollection>("GenVertex");
  writer.registerForWrite<fcc::GenJetCollection>("GenJet");
  // declaring new particles:
  TDatabasePDG* pdgBase = TDatabasePDG::Instance();

  double m=1.0;
  double g=1.0;
  pdgBase->AddParticle("N2","HNL", m, false, g, 0., "N2", 9900015);

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia8::Pythia pythia;
  pythia.readString("9900015:new = N2 N2 2 0 0 1.0 0.0 0.0 0.0 1.  0   1   0   1   0");



  pythia.readFile(card_file);
  //pythia.readString("23:oneChannel =  1 1 0 9900015 16");
  //pythia.SetParameters("9900015:oneChannel =  1 1 0 13 -13");



  //  std::cout<<"Listing Particle changes: "<<std::endl;
  //pythia.particleData.listAll();
  //std::cout<<"Done listing"<<std::endl;
  //std::cout<<"Now Reading new Table"<<std::endl;
  //  pythia.particleData.readXML("/afs/cern.ch/user/m/mchrzasz/FCC/fcc-physics/pythia8/ParticleData_HNL.xml");
  //std::cout<<"Changes in particles: "<<std::endl;
  //pythia.particleData.listChanged();

  //  pythia.particleData.listAll();

  //pythia.readFile("/afs/cern.ch/user/m/mchrzasz/FCC/fcc-physics/pythia8/ParticleData_HNL_short.xml");





  pythia.init();

  unsigned nevents = pythia.mode("Main:numberOfEvents");

  // Interface for conversion from Pythia8::Event to HepMC event.
  HepMC::Pythia8ToHepMC ToHepMC;
  Pythia8::ParticleData PDG;

  // Fastjet

   fastjet::Strategy strategy = fastjet::Best;
  fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;
  double rparam = 0.5;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, rparam,
                                 recomb_scheme, strategy);

  std::vector<fastjet::PseudoJet> input_particles;

  for(unsigned iev=0; iev<nevents; ++iev) {

    auto evinfo = evinfocoll.create();
    evinfo.Number(iev);

    if (!pythia.next()) continue;
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    hepmcevt->use_units(HepMC::Units::GEV, HepMC::Units::MM);
    ToHepMC.fill_next_event( pythia, hepmcevt );

    //    std::cout<<"Nvtx = "<<hepmcevt->vertices_size()<<std::endl;
    typedef std::map<HepMC::GenVertex*, fcc::GenVertex > VertexMap;
    VertexMap vtx_map;
    for ( HepMC::GenEvent::vertex_iterator iv = hepmcevt->vertices_begin();
	  iv != hepmcevt->vertices_end(); ++iv ) {
      const HepMC::FourVector& vpos = (*iv)->position();
      fcc::GenVertex vtx = vcoll.create();
      vtx.Position().X = vpos.x();
      vtx.Position().Y = vpos.y();
      vtx.Position().Z = vpos.z();
      vtx.Ctau(vpos.t());
      vtx_map.emplace(*iv, vtx);
    }
    input_particles.clear();
    for ( HepMC::GenEvent::particle_iterator ip = hepmcevt->particles_begin();
	  ip != hepmcevt->particles_end(); ++ip ) {
      HepMC::GenParticle* hepmcptc = *ip;
      fcc::MCParticle ptc = pcoll.create();
      fcc::BareParticle& core = ptc.Core();
      core.Type = hepmcptc->pdg_id();
      core.Charge = pythia.particleData.charge(core.Type);
      core.Status = hepmcptc->status();
      core.P4.Px = hepmcptc->momentum().px();
      core.P4.Py = hepmcptc->momentum().py();
      core.P4.Pz = hepmcptc->momentum().pz();
      core.P4.Mass = hepmcptc->momentum().m();

      if(core.Status==1) {
	input_particles.push_back( fastjet::PseudoJet(hepmcptc->momentum().px(),
						      hepmcptc->momentum().py(),
						      hepmcptc->momentum().pz(),
						      hepmcptc->momentum().e() ));
      }
      std::vector<fastjet::PseudoJet> input_particles;

      typedef VertexMap::const_iterator IVM;
      IVM prodvtx = vtx_map.find(hepmcptc->production_vertex());
      if(prodvtx!=vtx_map.end()) {
	ptc.StartVertex(prodvtx->second);
      }

      IVM endvtx = vtx_map.find(hepmcptc->end_vertex());
      if(endvtx!=vtx_map.end()) {
	ptc.EndVertex(endvtx->second);
      }

    }
    fastjet::ClusterSequence clust_seq(input_particles, jet_def);
    double ptmin = 10;
    std::vector<fastjet::PseudoJet> jets = clust_seq.inclusive_jets(ptmin);
    for(unsigned i=0; i<jets.size(); ++i) {
      fcc::GenJet genjet = genjetcoll.create();
      fcc::BareJet& core = genjet.Core();
      core.P4.Px = jets[i].px();
      core.P4.Py = jets[i].py();
      core.P4.Pz = jets[i].pz();
      core.P4.Mass = jets[i].m();
    }

    writer.writeEvent();
    store.clearCollections();
    delete hepmcevt;
  }

  writer.finish();

  return 0;
}


