/*

  CMS EXO-19-010 disappearing track search 139 fb^-1   
  -- arXiv:2004.05153  
  Recast By Mark Goodsell (goodsell@lpthe.jussieu.fr) 
  -- arXiv:2106.08815  

  CDFO Customization By Lucas Magno D. Ramos (lucasmdr@if.usp.br)
  -- arXiv:2404.16086
*/ 


#include "SampleAnalyzer/User/Analyzer/cms_exo_19_010.h"
using namespace MA5;
using namespace std;
												
MAbool LV_equal(MALorentzVector p, MALorentzVector q){
	if (p.X()==q.X() && p.Y()==q.Y() && p.Z()==q.Z() && p.T()==q.T()) return true;
	return false;
	}

												

void get_which_hits(const std::vector<double> &radii, std::vector<bool> &outhits, double maxR, double Rdec, double Rprod)
{
  //std::vector<bool> outhits;
  for(auto discR : radii)
  {
    if(discR > maxR) return;

    if(Rprod > discR)
    {
      outhits.push_back(false);
      continue;
    }

    if(Rdec < discR)
    {
      outhits.push_back(false);
      continue;
    }
    
    outhits.push_back(true);

  }

  
}

std::vector<bool> get_truth_hits(const RecTrackFormat* p)
{
  // This function to determine how many tracker layers are hit

  // strategy: Determine truth level and then apply per-hit efficiencies
  // Start with pixel detector and work outwards: TIB, TID, TEC, TOB
  // Pixel detector: https://lss.fnal.gov/archive/design/fermilab-design-2012-02.pdf figure 1.2, 2.1, table 2.1
  // Tracker: 1405.6569

  std::vector<bool> outhits;
  outhits.clear();

  
  //MAfloat64 abseta=fabs(p->etaCalo());
  MAfloat64 abseta=fabs(p->mc()->momentum().Eta());

  // don't bother with large eta
  if(abseta > 2.5)
  {
    //outhits.push_back(false);
    return outhits;
  }

  MALorentzVector vstart = p->ProductionVertex();

  double ppz = fabs(p->momentum().Pz());
  double ppT=p->momentum().Pt();
  
  double ootantheta=0.0;
  if(ppz > 0.0)
  {
    ootantheta = ppz/ppT; // nb for abseta > 2.5 we don't worry about zeroes here
  }

  double vdecz, vdecR;
  if(p->mc()->ctau() > 1.0e-3)
  {
    MALorentzVector vdec = p->mc()->decay_vertex();
    vdecz=fabs(vdec.Pz());
    vdecR=vdec.Pt();
  }
  else
  {
    vdecR=1.0e6;
    vdecz=vdecR*ootantheta;
  }

  // for now assume tracks start at z=0
  double vstartz=fabs(vstart.Pz());
  double vstartR=vstart.Pt();
  

  

  // Pixel
  // https://arxiv.org/pdf/0911.5434.pdf for 
  // Pixel barrel layers at e = 3, 6.8, 10.2, 16 cm
  // pixel z extends to |z| = 53.3/2 cm
  // This means |eta| < 1.3 is entirely inside the pixel barrel
  // Pre-upgrade: The endcap disks, extending from 6 to 15 cm in radius, are placed at z = ±35.5 cm and z = ±48.5 cm.
  // Post-upgrade: three endcap disks, see figure 2.1 in https://lss.fnal.gov/archive/design/fermilab-design-2012-02.pdf
  // let us guess that these are at 30, 40, 50 cm

  //1710.03842: Theradii of the four barrel layers are 2.9 cm, 6.8 cm, 10.9 cm, and 16.0 cm. 
  //The three forward disks arepositioned alongzat 3.2 cm, 3.9 cm, and 4.8 cm. 
  // Each disk is composed of two rings of moduleswith average radii of 12.8 cm and 7.8 cm.

  //static vector<double> pixelbarrelR= {30.0,68.0,102.0,160.0};
  static vector<double> pixelbarrelR= {29.0,68.0,109.0,160.0};  


  //static vector<double> pixelForwardZ = {300.0,400.0,500.0};
  static vector<double> pixelForwardZ = {320.0,390.0,480.0};


  // cannot escape inner barrel of pixel altogether
  double maxpixelR;
  if(ootantheta > 1)
  {
    maxpixelR=266.5/ootantheta;
  }
  else
  {
    maxpixelR=200.0;
  }

  get_which_hits(pixelbarrelR, outhits, maxpixelR, vdecR, vstartR);

  // check if we need pixel endcaps
  if(abseta > 1.3) 
  {
    
    double maxpixelZ=150.0*ootantheta;
    
    get_which_hits(pixelForwardZ, outhits, maxpixelZ, vdecz, vstartz);

  }
  
  if(outhits.size() < 4) 
  {
    //std::cout << "Track shorter than 4 pixel hits! " << abseta << std::endl;
    // Have I done something wrong for this case? Here I will pad the tracks with "false"
    // because we first check if there are any missing pixel hits
    for(int i=0; i < 4-outhits.size(); i++)
    {
      outhits.push_back(false);
    }

  }




  // Now for tracker barrel, see CMS 2008 report for some of these
  /*
  The Tracker Inner Barrel (TIB)
and Disks (TID) cover r < 55 cm and | z | < 118 cm, and are composed of four barrel layers,
supplemented by three disks at each end. The Tracker Outer Barrel (TOB) covers r > 55 cm and
| z | < 118 cm and consists of six barrel layers. The Tracker EndCaps (TEC) cover the region 124 <
| z | < 282 cm. Each TEC is composed of nine disks, each containing up to seven concentric
rings of silicon strip modules, yielding a range of resolutions similar to that of the TOB.
*/

  
  //static vector<double> TIBR= {230.0, 300.0, 400.0,500.0};
  // See CMS 2008 report page 64, these exetend -700 to 700 mm
  static vector<double> TIBR= {255.0, 339.0, 418.5,498.0};

  // The TID± are assemblies of three disks placed in z between ±800 mm and ±900 mm. The
  // disks are identical and each one consists of three rings which span the radius from roughly 200 mm
  // to 500 mm.
  // Here I use values inferred from Figure 3.34 in the CMS 2008 report
  static vector<double> TIDZ= {775.0,900.0,1025.0};

  // Close!!
  //static vector<double> TOBR = {600.0,700.0,800.0,890.0,980.0,1060.0};
  // See CMS 2008 report page 68
  static vector<double> TOBR = {608.0,692.0,780.0,868.0,965.0,1080.0};

  // Inferred from Figure 3.34 in the CMS 2008 report
  static vector<double> TECZ= {1250.0,1400.0,1550.0,1700.0,1950.0,2000.0,2225.0,2450.0,2700.0};

  // inner barrel
  double maxIBR;

  // inner barrel up to |z| = 580mm, r < 550
  //if(ootantheta < 1.05)
  // inner barrel up to |z| = 700mm, r < 550 -> 1/tanbeta < 700/550 = 1.27
  if(ootantheta < 1.27)
  {
    maxIBR=550.0;
    get_which_hits(TIBR, outhits, maxIBR, vdecR, vstartR);
  }
  else
  {
    maxIBR=550.0/ootantheta;
    get_which_hits(TIBR, outhits, maxIBR, vdecR, vstartR);
    //outer inner endcap, only possible if ootantheta > 1.05
    
    //double maxIBZ=580.0*ootantheta;
    double maxIBZ=700.0*ootantheta;
    get_which_hits(TIDZ, outhits, maxIBZ, vdecz, vstartz);
  }
  

  // now outer barrel, r > 55 cm and | z | < 118 cm 
  if(ootantheta < 2.145)
  {
    double maxOBR=1200.0;
    get_which_hits(TOBR, outhits, maxOBR, vdecR, vstartR);

  }
  else
  {
    double maxOBR=1180.0/ootantheta;
    get_which_hits(TOBR, outhits, maxOBR, vdecR, vstartR);

    // Now check the endcaps

    double maxTEZ=1180.0*ootantheta;
    get_which_hits(TECZ, outhits, maxTEZ, vdecz, vstartz);
  }



  return outhits;
}

class charged_track  {
   private:
    
   public:
    const RecTrackFormat* p; // in future change this to a vector or something else
    std::vector<bool> _hits;
    MAbool chdau_iso;
    
    charged_track() { }

    charged_track(const RecTrackFormat* q)
    {
      p = q;
      _hits = get_truth_hits(p);
      chdau_iso=true;
    };

    charged_track(const RecTrackFormat &q, std::mt19937 &engine, std::uniform_real_distribution<double> &rd)
    {
      p = &q;
      chdau_iso=true;
      std::vector<bool> newhits = get_truth_hits(p);

      // Now model efficiencies: first a per-hit efficiency
      // 94.5 percent is really for the earlier data; the later data claims 99% although I don't know whether to believe it
      int nhits=0;
      for(auto hit : newhits)
      {
        if(hit)
        {
          double prob = rd(engine);
          //cout << prob << ", " ;
          if(prob < 0.945)
          {
            _hits.push_back(true);
            nhits++;
            //std::cout << "1,";
            continue;
          }
          
        }
        _hits.push_back(false);
        //std::cout << "0,";
      }
      //cout << std::endl;

      // Next: we know that tracks with a smaller number of total hits 
      // see 1405.6569 page 21, "As a consequence, weaker selection criteria can be applied for tracks having many hit layers"
      // Apply 70% chance of reconstructing if nhits=4, 80% if nhits=5, 90% if 6 100% if > 6
      
      if(nhits < 7)
      {
        double recoprob=1.0;
        if(nhits == 4) recoprob=0.7;
        else if (nhits == 5) recoprob = 0.8;
        else if (nhits == 6) recoprob = 0.9;
        
        if(rd(engine) > recoprob)
        {
          _hits.clear();
          

        }

      }


   
    }

  // for some reason etaCalo is giving 0 
  //double abseta() { return fabs(p->etaCalo()); }
  double eta() {return p->mc()->momentum().Eta(); }
  double abseta() {return fabs(this->eta());}
  //double eta() {return p->etaCalo(); }
  double Pt() { return p->mc()->momentum().Pt();}

  double E() { return p->mc()->momentum().E(); }
  // This probably not initialised either
  //double phi() { return p->phiCalo(); }
  double phi() { return p->mc()->momentum().Phi(); }
  int abspid() {return abs(p->mc()->pdgid());}
  int pid() {return p->mc()->pdgid();}
  const MALorentzVector& momentum() const { return p->momentum(); }
  
};

																
//------------------------------------------------------------------------------------------------------------------------------

MAbool chdau_switch(std::vector<charged_track*> tracks){
	if (tracks.size()<1) return false;
	for (auto track : tracks){
		if (track->chdau_iso) return true;
	}
	return false;
}

//------------------------------------------------------------------------------------------------------------------------------
																


// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
bool cms_exo_19_010::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{
  cout << "BEGIN Initialization" << endl;

  PHYSICS->mcConfig().Reset();


    // definition of the multiparticle "hadronic"
    PHYSICS->mcConfig().AddHadronicId(-20543);
    PHYSICS->mcConfig().AddHadronicId(-20533);
    PHYSICS->mcConfig().AddHadronicId(-20523);
    PHYSICS->mcConfig().AddHadronicId(-20513);
    PHYSICS->mcConfig().AddHadronicId(-20433);
    PHYSICS->mcConfig().AddHadronicId(-20423);
    PHYSICS->mcConfig().AddHadronicId(-20413);
    PHYSICS->mcConfig().AddHadronicId(-20323);
    PHYSICS->mcConfig().AddHadronicId(-20313);
    PHYSICS->mcConfig().AddHadronicId(-20213);
    PHYSICS->mcConfig().AddHadronicId(-10543);
    PHYSICS->mcConfig().AddHadronicId(-10541);
    PHYSICS->mcConfig().AddHadronicId(-10533);
    PHYSICS->mcConfig().AddHadronicId(-10531);
    PHYSICS->mcConfig().AddHadronicId(-10523);
    PHYSICS->mcConfig().AddHadronicId(-10521);
    PHYSICS->mcConfig().AddHadronicId(-10513);
    PHYSICS->mcConfig().AddHadronicId(-10511);
    PHYSICS->mcConfig().AddHadronicId(-10433);
    PHYSICS->mcConfig().AddHadronicId(-10431);
    PHYSICS->mcConfig().AddHadronicId(-10423);
    PHYSICS->mcConfig().AddHadronicId(-10421);
    PHYSICS->mcConfig().AddHadronicId(-10413);
    PHYSICS->mcConfig().AddHadronicId(-10411);
    PHYSICS->mcConfig().AddHadronicId(-10323);
    PHYSICS->mcConfig().AddHadronicId(-10321);
    PHYSICS->mcConfig().AddHadronicId(-10313);
    PHYSICS->mcConfig().AddHadronicId(-10311);
    PHYSICS->mcConfig().AddHadronicId(-10213);
    PHYSICS->mcConfig().AddHadronicId(-10211);
    PHYSICS->mcConfig().AddHadronicId(-5554);
    PHYSICS->mcConfig().AddHadronicId(-5544);
    PHYSICS->mcConfig().AddHadronicId(-5542);
    PHYSICS->mcConfig().AddHadronicId(-5534);
    PHYSICS->mcConfig().AddHadronicId(-5532);
    PHYSICS->mcConfig().AddHadronicId(-5524);
    PHYSICS->mcConfig().AddHadronicId(-5522);
    PHYSICS->mcConfig().AddHadronicId(-5514);
    PHYSICS->mcConfig().AddHadronicId(-5512);
    PHYSICS->mcConfig().AddHadronicId(-5503);
    PHYSICS->mcConfig().AddHadronicId(-5444);
    PHYSICS->mcConfig().AddHadronicId(-5442);
    PHYSICS->mcConfig().AddHadronicId(-5434);
    PHYSICS->mcConfig().AddHadronicId(-5432);
    PHYSICS->mcConfig().AddHadronicId(-5424);
    PHYSICS->mcConfig().AddHadronicId(-5422);
    PHYSICS->mcConfig().AddHadronicId(-5414);
    PHYSICS->mcConfig().AddHadronicId(-5412);
    PHYSICS->mcConfig().AddHadronicId(-5403);
    PHYSICS->mcConfig().AddHadronicId(-5401);
    PHYSICS->mcConfig().AddHadronicId(-5342);
    PHYSICS->mcConfig().AddHadronicId(-5334);
    PHYSICS->mcConfig().AddHadronicId(-5332);
    PHYSICS->mcConfig().AddHadronicId(-5324);
    PHYSICS->mcConfig().AddHadronicId(-5322);
    PHYSICS->mcConfig().AddHadronicId(-5314);
    PHYSICS->mcConfig().AddHadronicId(-5312);
    PHYSICS->mcConfig().AddHadronicId(-5303);
    PHYSICS->mcConfig().AddHadronicId(-5301);
    PHYSICS->mcConfig().AddHadronicId(-5242);
    PHYSICS->mcConfig().AddHadronicId(-5232);
    PHYSICS->mcConfig().AddHadronicId(-5224);
    PHYSICS->mcConfig().AddHadronicId(-5222);
    PHYSICS->mcConfig().AddHadronicId(-5214);
    PHYSICS->mcConfig().AddHadronicId(-5212);
    PHYSICS->mcConfig().AddHadronicId(-5203);
    PHYSICS->mcConfig().AddHadronicId(-5201);
    PHYSICS->mcConfig().AddHadronicId(-5142);
    PHYSICS->mcConfig().AddHadronicId(-5132);
    PHYSICS->mcConfig().AddHadronicId(-5122);
    PHYSICS->mcConfig().AddHadronicId(-5114);
    PHYSICS->mcConfig().AddHadronicId(-5112);
    PHYSICS->mcConfig().AddHadronicId(-5103);
    PHYSICS->mcConfig().AddHadronicId(-5101);
    PHYSICS->mcConfig().AddHadronicId(-4444);
    PHYSICS->mcConfig().AddHadronicId(-4434);
    PHYSICS->mcConfig().AddHadronicId(-4432);
    PHYSICS->mcConfig().AddHadronicId(-4424);
    PHYSICS->mcConfig().AddHadronicId(-4422);
    PHYSICS->mcConfig().AddHadronicId(-4414);
    PHYSICS->mcConfig().AddHadronicId(-4412);
    PHYSICS->mcConfig().AddHadronicId(-4403);
    PHYSICS->mcConfig().AddHadronicId(-4334);
    PHYSICS->mcConfig().AddHadronicId(-4332);
    PHYSICS->mcConfig().AddHadronicId(-4324);
    PHYSICS->mcConfig().AddHadronicId(-4322);
    PHYSICS->mcConfig().AddHadronicId(-4314);
    PHYSICS->mcConfig().AddHadronicId(-4312);
    PHYSICS->mcConfig().AddHadronicId(-4303);
    PHYSICS->mcConfig().AddHadronicId(-4301);
    PHYSICS->mcConfig().AddHadronicId(-4232);
    PHYSICS->mcConfig().AddHadronicId(-4224);
    PHYSICS->mcConfig().AddHadronicId(-4222);
    PHYSICS->mcConfig().AddHadronicId(-4214);
    PHYSICS->mcConfig().AddHadronicId(-4212);
    PHYSICS->mcConfig().AddHadronicId(-4203);
    PHYSICS->mcConfig().AddHadronicId(-4201);
    PHYSICS->mcConfig().AddHadronicId(-4132);
    PHYSICS->mcConfig().AddHadronicId(-4122);
    PHYSICS->mcConfig().AddHadronicId(-4114);
    PHYSICS->mcConfig().AddHadronicId(-4112);
    PHYSICS->mcConfig().AddHadronicId(-4103);
    PHYSICS->mcConfig().AddHadronicId(-4101);
    PHYSICS->mcConfig().AddHadronicId(-3334);
    PHYSICS->mcConfig().AddHadronicId(-3324);
    PHYSICS->mcConfig().AddHadronicId(-3322);
    PHYSICS->mcConfig().AddHadronicId(-3314);
    PHYSICS->mcConfig().AddHadronicId(-3312);
    PHYSICS->mcConfig().AddHadronicId(-3303);
    PHYSICS->mcConfig().AddHadronicId(-3224);
    PHYSICS->mcConfig().AddHadronicId(-3222);
    PHYSICS->mcConfig().AddHadronicId(-3214);
    PHYSICS->mcConfig().AddHadronicId(-3212);
    PHYSICS->mcConfig().AddHadronicId(-3203);
    PHYSICS->mcConfig().AddHadronicId(-3201);
    PHYSICS->mcConfig().AddHadronicId(-3122);
    PHYSICS->mcConfig().AddHadronicId(-3114);
    PHYSICS->mcConfig().AddHadronicId(-3112);
    PHYSICS->mcConfig().AddHadronicId(-3103);
    PHYSICS->mcConfig().AddHadronicId(-3101);
    PHYSICS->mcConfig().AddHadronicId(-2224);
    PHYSICS->mcConfig().AddHadronicId(-2214);
    PHYSICS->mcConfig().AddHadronicId(-2212);
    PHYSICS->mcConfig().AddHadronicId(-2203);
    PHYSICS->mcConfig().AddHadronicId(-2114);
    PHYSICS->mcConfig().AddHadronicId(-2112);
    PHYSICS->mcConfig().AddHadronicId(-2103);
    PHYSICS->mcConfig().AddHadronicId(-2101);
    PHYSICS->mcConfig().AddHadronicId(-1114);
    PHYSICS->mcConfig().AddHadronicId(-1103);
    PHYSICS->mcConfig().AddHadronicId(-545);
    PHYSICS->mcConfig().AddHadronicId(-543);
    PHYSICS->mcConfig().AddHadronicId(-541);
    PHYSICS->mcConfig().AddHadronicId(-535);
    PHYSICS->mcConfig().AddHadronicId(-533);
    PHYSICS->mcConfig().AddHadronicId(-531);
    PHYSICS->mcConfig().AddHadronicId(-525);
    PHYSICS->mcConfig().AddHadronicId(-523);
    PHYSICS->mcConfig().AddHadronicId(-521);
    PHYSICS->mcConfig().AddHadronicId(-515);
    PHYSICS->mcConfig().AddHadronicId(-513);
    PHYSICS->mcConfig().AddHadronicId(-511);
    PHYSICS->mcConfig().AddHadronicId(-435);
    PHYSICS->mcConfig().AddHadronicId(-433);
    PHYSICS->mcConfig().AddHadronicId(-431);
    PHYSICS->mcConfig().AddHadronicId(-425);
    PHYSICS->mcConfig().AddHadronicId(-423);
    PHYSICS->mcConfig().AddHadronicId(-421);
    PHYSICS->mcConfig().AddHadronicId(-415);
    PHYSICS->mcConfig().AddHadronicId(-413);
    PHYSICS->mcConfig().AddHadronicId(-411);
    PHYSICS->mcConfig().AddHadronicId(-325);
    PHYSICS->mcConfig().AddHadronicId(-323);
    PHYSICS->mcConfig().AddHadronicId(-321);
    PHYSICS->mcConfig().AddHadronicId(-315);
    PHYSICS->mcConfig().AddHadronicId(-313);
    PHYSICS->mcConfig().AddHadronicId(-311);
    PHYSICS->mcConfig().AddHadronicId(-215);
    PHYSICS->mcConfig().AddHadronicId(-213);
    PHYSICS->mcConfig().AddHadronicId(-211);
    PHYSICS->mcConfig().AddHadronicId(111);
    PHYSICS->mcConfig().AddHadronicId(113);
    PHYSICS->mcConfig().AddHadronicId(115);
    PHYSICS->mcConfig().AddHadronicId(130);
    PHYSICS->mcConfig().AddHadronicId(211);
    PHYSICS->mcConfig().AddHadronicId(213);
    PHYSICS->mcConfig().AddHadronicId(215);
    PHYSICS->mcConfig().AddHadronicId(221);
    PHYSICS->mcConfig().AddHadronicId(223);
    PHYSICS->mcConfig().AddHadronicId(225);
    PHYSICS->mcConfig().AddHadronicId(310);
    PHYSICS->mcConfig().AddHadronicId(311);
    PHYSICS->mcConfig().AddHadronicId(313);
    PHYSICS->mcConfig().AddHadronicId(315);
    PHYSICS->mcConfig().AddHadronicId(321);
    PHYSICS->mcConfig().AddHadronicId(323);
    PHYSICS->mcConfig().AddHadronicId(325);
    PHYSICS->mcConfig().AddHadronicId(331);
    PHYSICS->mcConfig().AddHadronicId(333);
    PHYSICS->mcConfig().AddHadronicId(335);
    PHYSICS->mcConfig().AddHadronicId(411);
    PHYSICS->mcConfig().AddHadronicId(413);
    PHYSICS->mcConfig().AddHadronicId(415);
    PHYSICS->mcConfig().AddHadronicId(421);
    PHYSICS->mcConfig().AddHadronicId(423);
    PHYSICS->mcConfig().AddHadronicId(425);
    PHYSICS->mcConfig().AddHadronicId(431);
    PHYSICS->mcConfig().AddHadronicId(433);
    PHYSICS->mcConfig().AddHadronicId(435);
    PHYSICS->mcConfig().AddHadronicId(441);
    PHYSICS->mcConfig().AddHadronicId(443);
    PHYSICS->mcConfig().AddHadronicId(445);
    PHYSICS->mcConfig().AddHadronicId(511);
    PHYSICS->mcConfig().AddHadronicId(513);
    PHYSICS->mcConfig().AddHadronicId(515);
    PHYSICS->mcConfig().AddHadronicId(521);
    PHYSICS->mcConfig().AddHadronicId(523);
    PHYSICS->mcConfig().AddHadronicId(525);
    PHYSICS->mcConfig().AddHadronicId(531);
    PHYSICS->mcConfig().AddHadronicId(533);
    PHYSICS->mcConfig().AddHadronicId(535);
    PHYSICS->mcConfig().AddHadronicId(541);
    PHYSICS->mcConfig().AddHadronicId(543);
    PHYSICS->mcConfig().AddHadronicId(545);
    PHYSICS->mcConfig().AddHadronicId(551);
    PHYSICS->mcConfig().AddHadronicId(553);
    PHYSICS->mcConfig().AddHadronicId(555);
    PHYSICS->mcConfig().AddHadronicId(1103);
    PHYSICS->mcConfig().AddHadronicId(1114);
    PHYSICS->mcConfig().AddHadronicId(2101);
    PHYSICS->mcConfig().AddHadronicId(2103);
    PHYSICS->mcConfig().AddHadronicId(2112);
    PHYSICS->mcConfig().AddHadronicId(2114);
    PHYSICS->mcConfig().AddHadronicId(2203);
    PHYSICS->mcConfig().AddHadronicId(2212);
    PHYSICS->mcConfig().AddHadronicId(2214);
    PHYSICS->mcConfig().AddHadronicId(2224);
    PHYSICS->mcConfig().AddHadronicId(3101);
    PHYSICS->mcConfig().AddHadronicId(3103);
    PHYSICS->mcConfig().AddHadronicId(3112);
    PHYSICS->mcConfig().AddHadronicId(3114);
    PHYSICS->mcConfig().AddHadronicId(3122);
    PHYSICS->mcConfig().AddHadronicId(3201);
    PHYSICS->mcConfig().AddHadronicId(3203);
    PHYSICS->mcConfig().AddHadronicId(3212);
    PHYSICS->mcConfig().AddHadronicId(3214);
    PHYSICS->mcConfig().AddHadronicId(3222);
    PHYSICS->mcConfig().AddHadronicId(3224);
    PHYSICS->mcConfig().AddHadronicId(3303);
    PHYSICS->mcConfig().AddHadronicId(3312);
    PHYSICS->mcConfig().AddHadronicId(3314);
    PHYSICS->mcConfig().AddHadronicId(3322);
    PHYSICS->mcConfig().AddHadronicId(3324);
    PHYSICS->mcConfig().AddHadronicId(3334);
    PHYSICS->mcConfig().AddHadronicId(4101);
    PHYSICS->mcConfig().AddHadronicId(4103);
    PHYSICS->mcConfig().AddHadronicId(4112);
    PHYSICS->mcConfig().AddHadronicId(4114);
    PHYSICS->mcConfig().AddHadronicId(4122);
    PHYSICS->mcConfig().AddHadronicId(4132);
    PHYSICS->mcConfig().AddHadronicId(4201);
    PHYSICS->mcConfig().AddHadronicId(4203);
    PHYSICS->mcConfig().AddHadronicId(4212);
    PHYSICS->mcConfig().AddHadronicId(4214);
    PHYSICS->mcConfig().AddHadronicId(4222);
    PHYSICS->mcConfig().AddHadronicId(4224);
    PHYSICS->mcConfig().AddHadronicId(4232);
    PHYSICS->mcConfig().AddHadronicId(4301);
    PHYSICS->mcConfig().AddHadronicId(4303);
    PHYSICS->mcConfig().AddHadronicId(4312);
    PHYSICS->mcConfig().AddHadronicId(4314);
    PHYSICS->mcConfig().AddHadronicId(4322);
    PHYSICS->mcConfig().AddHadronicId(4324);
    PHYSICS->mcConfig().AddHadronicId(4332);
    PHYSICS->mcConfig().AddHadronicId(4334);
    PHYSICS->mcConfig().AddHadronicId(4403);
    PHYSICS->mcConfig().AddHadronicId(4412);
    PHYSICS->mcConfig().AddHadronicId(4414);
    PHYSICS->mcConfig().AddHadronicId(4422);
    PHYSICS->mcConfig().AddHadronicId(4424);
    PHYSICS->mcConfig().AddHadronicId(4432);
    PHYSICS->mcConfig().AddHadronicId(4434);
    PHYSICS->mcConfig().AddHadronicId(4444);
    PHYSICS->mcConfig().AddHadronicId(5101);
    PHYSICS->mcConfig().AddHadronicId(5103);
    PHYSICS->mcConfig().AddHadronicId(5112);
    PHYSICS->mcConfig().AddHadronicId(5114);
    PHYSICS->mcConfig().AddHadronicId(5122);
    PHYSICS->mcConfig().AddHadronicId(5132);
    PHYSICS->mcConfig().AddHadronicId(5142);
    PHYSICS->mcConfig().AddHadronicId(5201);
    PHYSICS->mcConfig().AddHadronicId(5203);
    PHYSICS->mcConfig().AddHadronicId(5212);
    PHYSICS->mcConfig().AddHadronicId(5214);
    PHYSICS->mcConfig().AddHadronicId(5222);
    PHYSICS->mcConfig().AddHadronicId(5224);
    PHYSICS->mcConfig().AddHadronicId(5232);
    PHYSICS->mcConfig().AddHadronicId(5242);
    PHYSICS->mcConfig().AddHadronicId(5301);
    PHYSICS->mcConfig().AddHadronicId(5303);
    PHYSICS->mcConfig().AddHadronicId(5312);
    PHYSICS->mcConfig().AddHadronicId(5314);
    PHYSICS->mcConfig().AddHadronicId(5322);
    PHYSICS->mcConfig().AddHadronicId(5324);
    PHYSICS->mcConfig().AddHadronicId(5332);
    PHYSICS->mcConfig().AddHadronicId(5334);
    PHYSICS->mcConfig().AddHadronicId(5342);
    PHYSICS->mcConfig().AddHadronicId(5401);
    PHYSICS->mcConfig().AddHadronicId(5403);
    PHYSICS->mcConfig().AddHadronicId(5412);
    PHYSICS->mcConfig().AddHadronicId(5414);
    PHYSICS->mcConfig().AddHadronicId(5422);
    PHYSICS->mcConfig().AddHadronicId(5424);
    PHYSICS->mcConfig().AddHadronicId(5432);
    PHYSICS->mcConfig().AddHadronicId(5434);
    PHYSICS->mcConfig().AddHadronicId(5442);
    PHYSICS->mcConfig().AddHadronicId(5444);
    PHYSICS->mcConfig().AddHadronicId(5503);
    PHYSICS->mcConfig().AddHadronicId(5512);
    PHYSICS->mcConfig().AddHadronicId(5514);
    PHYSICS->mcConfig().AddHadronicId(5522);
    PHYSICS->mcConfig().AddHadronicId(5524);
    PHYSICS->mcConfig().AddHadronicId(5532);
    PHYSICS->mcConfig().AddHadronicId(5534);
    PHYSICS->mcConfig().AddHadronicId(5542);
    PHYSICS->mcConfig().AddHadronicId(5544);
    PHYSICS->mcConfig().AddHadronicId(5554);
    PHYSICS->mcConfig().AddHadronicId(10111);
    PHYSICS->mcConfig().AddHadronicId(10113);
    PHYSICS->mcConfig().AddHadronicId(10211);
    PHYSICS->mcConfig().AddHadronicId(10213);
    PHYSICS->mcConfig().AddHadronicId(10221);
    PHYSICS->mcConfig().AddHadronicId(10223);
    PHYSICS->mcConfig().AddHadronicId(10311);
    PHYSICS->mcConfig().AddHadronicId(10313);
    PHYSICS->mcConfig().AddHadronicId(10321);
    PHYSICS->mcConfig().AddHadronicId(10323);
    PHYSICS->mcConfig().AddHadronicId(10331);
    PHYSICS->mcConfig().AddHadronicId(10333);
    PHYSICS->mcConfig().AddHadronicId(10411);
    PHYSICS->mcConfig().AddHadronicId(10413);
    PHYSICS->mcConfig().AddHadronicId(10421);
    PHYSICS->mcConfig().AddHadronicId(10423);
    PHYSICS->mcConfig().AddHadronicId(10431);
    PHYSICS->mcConfig().AddHadronicId(10433);
    PHYSICS->mcConfig().AddHadronicId(10441);
    PHYSICS->mcConfig().AddHadronicId(10443);
    PHYSICS->mcConfig().AddHadronicId(10511);
    PHYSICS->mcConfig().AddHadronicId(10513);
    PHYSICS->mcConfig().AddHadronicId(10521);
    PHYSICS->mcConfig().AddHadronicId(10523);
    PHYSICS->mcConfig().AddHadronicId(10531);
    PHYSICS->mcConfig().AddHadronicId(10533);
    PHYSICS->mcConfig().AddHadronicId(10541);
    PHYSICS->mcConfig().AddHadronicId(10543);
    PHYSICS->mcConfig().AddHadronicId(10551);
    PHYSICS->mcConfig().AddHadronicId(10553);
    PHYSICS->mcConfig().AddHadronicId(20113);
    PHYSICS->mcConfig().AddHadronicId(20213);
    PHYSICS->mcConfig().AddHadronicId(20223);
    PHYSICS->mcConfig().AddHadronicId(20313);
    PHYSICS->mcConfig().AddHadronicId(20323);
    PHYSICS->mcConfig().AddHadronicId(20333);
    PHYSICS->mcConfig().AddHadronicId(20413);
    PHYSICS->mcConfig().AddHadronicId(20423);
    PHYSICS->mcConfig().AddHadronicId(20433);
    PHYSICS->mcConfig().AddHadronicId(20443);
    PHYSICS->mcConfig().AddHadronicId(20513);
    PHYSICS->mcConfig().AddHadronicId(20523);
    PHYSICS->mcConfig().AddHadronicId(20533);
    PHYSICS->mcConfig().AddHadronicId(20543);
    PHYSICS->mcConfig().AddHadronicId(20553);
    PHYSICS->mcConfig().AddHadronicId(100443);
    PHYSICS->mcConfig().AddHadronicId(100553);
    PHYSICS->mcConfig().AddHadronicId(9900440);
    PHYSICS->mcConfig().AddHadronicId(9900441);
    PHYSICS->mcConfig().AddHadronicId(9900443);
    PHYSICS->mcConfig().AddHadronicId(9900551);
    PHYSICS->mcConfig().AddHadronicId(9900553);
    PHYSICS->mcConfig().AddHadronicId(9910441);
    PHYSICS->mcConfig().AddHadronicId(9910551);

    // definition of the multiparticle "invisible"
    PHYSICS->mcConfig().AddInvisibleId(-16);
    PHYSICS->mcConfig().AddInvisibleId(-14);
    PHYSICS->mcConfig().AddInvisibleId(-12);
    PHYSICS->mcConfig().AddInvisibleId(12);
    PHYSICS->mcConfig().AddInvisibleId(14);
    PHYSICS->mcConfig().AddInvisibleId(16);
    PHYSICS->mcConfig().AddInvisibleId(1000022);
    PHYSICS->mcConfig().AddInvisibleId(1000039);

    // Initializing PhysicsService for RECO
    PHYSICS->recConfig().Reset();
  
  // initialize variables, histos

  cout << "--------------------------------------------------" <<std::endl;
  cout << "-- CMS disappearing track search 139 fb^-1      --" <<std::endl;
  cout << "-- arXiv:2004.05153                             --" << std::endl;
  cout << "-- By Mark Goodsell (goodsell@lpthe.jussieu.fr) --" << std::endl;
  cout << "--------------------------------------------------" <<std::endl;
  this->engine = std::mt19937(time(NULL));

  
  const std::vector<std::string> vecallSRs={"SR3_2015","SR3_2016A","SR3_2016B","SR1_2017","SR1_2018A","SR1_2018B","SR2_2017","SR2_2018A","SR2_2018B","SR3_2017","SR3_2018A","SR3_2018B",
  						"SR3_2015_dR","SR3_2016A_dR","SR3_2016B_dR","SR1_2017_dR","SR1_2018A_dR","SR1_2018B_dR","SR2_2017_dR","SR2_2018A_dR","SR2_2018B_dR","SR3_2017_dR","SR3_2018A_dR","SR3_2018B_dR"};

    std::string allSRs0[]={"SR3_2015","SR3_2016A","SR3_2016B","SR1_2017","SR1_2018A","SR1_2018B","SR2_2017","SR2_2018A","SR2_2018B","SR3_2017","SR3_2018A","SR3_2018B"};

    std::string allSRs[]={"SR3_2015","SR3_2016A","SR3_2016B","SR1_2017","SR1_2018A","SR1_2018B","SR2_2017","SR2_2018A","SR2_2018B","SR3_2017","SR3_2018A","SR3_2018B",
    "SR3_2015_dR","SR3_2016A_dR","SR3_2016B_dR","SR1_2017_dR","SR1_2018A_dR","SR1_2018B_dR","SR2_2017_dR","SR2_2018A_dR","SR2_2018B_dR","SR3_2017_dR","SR3_2018A_dR","SR3_2018B_dR"};

    std::string allSRs2[]={"SR1_2017","SR1_2018A","SR1_2018B","SR2_2017","SR2_2018A","SR2_2018B","SR3_2017","SR3_2018A","SR3_2018B",
    "SR1_2017_dR","SR1_2018A_dR","SR1_2018B_dR","SR2_2017_dR","SR2_2018A_dR","SR2_2018B_dR","SR3_2017_dR","SR3_2018A_dR","SR3_2018B_dR"};

    std::string alldRs[]={"SR3_2015_dR","SR3_2016A_dR","SR3_2016B_dR","SR1_2017_dR","SR1_2018A_dR","SR1_2018B_dR","SR2_2017_dR","SR2_2018A_dR","SR2_2018B_dR","SR3_2017_dR","SR3_2018A_dR","SR3_2018B_dR"};

    std::string all2015[]={"SR3_2015","SR3_2015_dR"};
    
    std::string all2016[]={"SR3_2016A","SR3_2016B"};
    std::string all2016dR[]={"SR3_2016A_dR","SR3_2016B_dR"};
   
    std::string all2017[]={"SR1_2017","SR2_2017","SR3_2017"};
    std::string all2017dR[]={"SR1_2017_dR","SR2_2017_dR","SR3_2017_dR"};
   
    std::string all2018A[]={"SR1_2018A","SR2_2018A","SR3_2018A"};
    std::string all2018AdR[]={"SR1_2018A_dR","SR2_2018A_dR","SR3_2018A_dR"};
   
    std::string all2018B[]={"SR1_2018B","SR2_2018B","SR3_2018B"};
    std::string all2018BdR[]={"SR1_2018B_dR","SR2_2018B_dR","SR3_2018B_dR"};




    for(std::string region : vecallSRs)
      {
	Manager()->AddRegionSelection(region);
      }
  

 

    Manager()->AddCut("MET2015",all2015);
    Manager()->AddCut("MET",allSRs2);

    Manager()->AddCut("OneGoodJet",allSRs);

    Manager()->AddCut("JetAngleCuts",allSRs);
    Manager()->AddCut("|Delta phi(leading jet, pTmiss)| > 0.5",allSRs);

    Manager()->AddCut(">=1 track with |eta| < 2.1",allSRs);
    Manager()->AddCut(">=1 track with p_T > 55",allSRs);
    Manager()->AddCut(">=1 track passing fiducial selections",allSRs);
    Manager()->AddCut(">=1 track with relative isolation",allSRs);
    Manager()->AddCut(">=1 track with > 3 or 4 pixel hits",allSRs);
    Manager()->AddCut(">=1 track with no missing inner/middle hits",allSRs);
    Manager()->AddCut(">=1 track with d0 < 0.2 mm",allSRs);
    Manager()->AddCut(">=1 track with dz < 5.0 mm",allSRs);

    Manager()->AddCut(">=1 track with DR(track,daughter)>0.2",alldRs);

    Manager()->AddCut(">=1 track with DR(track, jet) > 0.5",allSRs0);
    Manager()->AddCut(">=1 track with DR(track, jet) > 0.5.",alldRs);
    Manager()->AddCut(">=1 track with DR(track, electron) > 0.15",allSRs0);
    Manager()->AddCut(">=1 track with DR(track, electron) > 0.15.",alldRs);
    Manager()->AddCut(">=1 track with DR(track, muon) > 0.15",allSRs0);
    Manager()->AddCut(">=1 track with DR(track, muon) > 0.15.",alldRs);
    Manager()->AddCut(">=1 track with DR(track, tau) > 0.15",allSRs0);
    Manager()->AddCut(">=1 track with DR(track, tau) > 0.15.",alldRs);
    Manager()->AddCut(">=1 track with Ecalo < 10",allSRs0);
    Manager()->AddCut(">=1 track with Ecalo < 10.",alldRs);
    Manager()->AddCut(">=1 track with >=3 missing outer hits",allSRs0);
    Manager()->AddCut(">=1 track with >=3 missing outer hits.",alldRs);

    Manager()->AddCut("HEMVeto",all2018B);
    Manager()->AddCut("HEMVeto.",all2018BdR);

    Manager()->AddCut(">5 Layers 2015","SR3_2015");
    Manager()->AddCut(">5 Layers 2015.","SR3_2015_dR");
    Manager()->AddCut(">5 Layers 2016",all2016);
    Manager()->AddCut(">5 Layers 2016.",all2016dR);

    Manager()->AddCut("4 Layers 2017","SR1_2017");
    Manager()->AddCut("4 Layers 2017.","SR1_2017_dR");

    Manager()->AddCut("4 Layers 2018A","SR1_2018A");
    Manager()->AddCut("4 Layers 2018A.","SR1_2018A_dR");

    Manager()->AddCut("4 Layers 2018B","SR1_2018B");
    Manager()->AddCut("4 Layers 2018B.","SR1_2018B_dR");

    Manager()->AddCut("5 Layers 2017","SR2_2017");
    Manager()->AddCut("5 Layers 2017.","SR2_2017_dR");
    Manager()->AddCut("5 Layers 2018A","SR2_2018A");
    Manager()->AddCut("5 Layers 2018A.","SR2_2018A_dR");
    Manager()->AddCut("5 Layers 2018B","SR2_2018B");
    Manager()->AddCut("5 Layers 2018B.","SR2_2018B_dR");

    Manager()->AddCut(">5 Layers 2017","SR3_2017");
    Manager()->AddCut(">5 Layers 2017.","SR3_2017_dR");
    Manager()->AddCut(">5 Layers 2018A","SR3_2018A");
    Manager()->AddCut(">5 Layers 2018A.","SR3_2018A_dR");
    Manager()->AddCut(">5 Layers 2018B","SR3_2018B");
    Manager()->AddCut(">5 Layers 2018B.","SR3_2018B_dR");

    Manager()->AddCut("Reweight 2015","SR3_2015");
    Manager()->AddCut("Reweight 2015.","SR3_2015_dR");
    Manager()->AddCut("Reweight 2016A","SR3_2016A");
    Manager()->AddCut("Reweight 2016A.","SR3_2016A_dR");
    Manager()->AddCut("Reweight 2016B","SR3_2016B");
    Manager()->AddCut("Reweight 2016B.","SR3_2016B_dR");
    Manager()->AddCut("Reweight 2017",all2017);
    Manager()->AddCut("Reweight 2017.",all2017dR);
    Manager()->AddCut("Reweight 2018A",all2018A);
    Manager()->AddCut("Reweight 2018A.",all2018AdR);
    Manager()->AddCut("Reweight 2018B",all2018B);
    Manager()->AddCut("Reweight 2018B.",all2018BdR);

    double totallumi=2.7+8.3+27.4+42.0+21.0+39.0;
    
    lumiratios["2015"] = 2.7/totallumi;
    lumiratios["2016A"] = 8.3/totallumi;
    lumiratios["2016B"] = 27.4/totallumi;
    lumiratios["2017"] = 42.0/totallumi;
    lumiratios["2018A"] = 21.0/totallumi;
    lumiratios["2018B"] = 39.0/totallumi;
  
  cout << "END   Initialization" << endl;
  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void cms_exo_19_010::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
  cout << "BEGIN Finalization" << endl;
  // saving histos
  cout << "END   Finalization" << endl;
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool cms_exo_19_010::Execute(SampleFormat& sample, const EventFormat& event)
{

  
  double myWeight=1.;
  if (!Configuration().IsNoEventWeight()) myWeight=event.mc()->weight();
  else if(event.mc()->weight()!=0.) myWeight=event.mc()->weight();
  else
  {
    //////WARNING << "Found one event with a zero weight. Skipping...\n";
    return false;
  }
  
  
  Manager()->InitializeForNewEvent(myWeight);

  // The event loop starts here
  if(event.rec() ==0) return false;
  
  MALorentzVector pTmiss,pTmissNoMu;

  std::vector<const RecTrackFormat*> charginos;

  // Count signal electrons, muons, jets etc 

static std::uniform_real_distribution<double> rd(0.0,1.0);
bool chargino50=false;

std::vector<charged_track*> tracks;

MALorentzVector muonp;
MALorentzVector charginop;
 muonp.clear();
 charginop.clear();
 
 for(auto muon : event.rec()->muons())
  {
    muonp+=muon.momentum();

  }

double maxtrackpt=0.0;
double chpt=0.0;


 for(auto &track : event.rec()->tracks()) 
  {
    
    const MCParticleFormat* chargino = track.mc();

    // NOT for tracks that are final states: MA5 will have already included them in the ET calculation
    //if ( (PHYSICS->Id->IsFinalState(chargino)) || (chargino->decay_vertex().Pt() > 7.0e3) || (fabs(chargino->decay_vertex().Pz()) > 1.1e4))
    if ( (chargino->decay_vertex().Pt() > 7.0e3) || (fabs(chargino->decay_vertex().Pz()) > 1.1e4))
    {
      if(!(PHYSICS->Id->IsHadronic(chargino)))
	{
      // gets reconstructed as a muon 
	  charginop+=chargino->momentum(); 
      //new_fake_muons.push_back(chargino);
	}
    }
    
    chpt=chargino->momentum().Pt();

    if(chpt > maxtrackpt) maxtrackpt=chpt; 
     
    //if(fabs(track.etaCalo()) < 2.5)
    if(fabs(chargino->momentum().Eta()) < 2.5)
            { 
              charged_track* newtrack = new charged_track(track,engine,rd); 
              tracks.push_back(newtrack); 

              if((!chargino50) && (chargino->momentum().Pt() >50.0)) chargino50=true; 
               
 
            } 
 


  }

 //////// NOW DO MET CALC! 

 //pTmiss = event.rec()->MET().momentum();
 // pTmissNoMu=pTmiss+muonp;

 pTmissNoMu = event.rec()->MET().momentum();
 
 pTmiss=pTmissNoMu-charginop;

 pTmissNoMu =  pTmissNoMu +muonp;
 double MET=pTmiss.Pt();
 double METnomu=pTmissNoMu.Pt();



 bool passMET=false; 
  bool passMET2015=false; 

 if((chargino50) && (MET > 105.0)) passMET=true; 
  if(MET > 120.0) passMET = true; 
  if(METnomu < 120.0)  
  { 
    passMET = false; 
  } 
  else 
  { 
    passMET = true; 
  } 
 
  if((chargino50) && (MET > 75.0) ) passMET2015=true; 
  if((MET > 90.0) || (METnomu > 90.0)) passMET2015 = true; 
  if(MET < 100.0) passMET2015 = false; 
 
   
  double eff=1.0; 
  //  1903.06078 figure 5, since the first cut on MET is the trigger. It uses HLT with 90 GeV or 120 GeV 
  if(passMET2015) 
  { 
 
    double tMET=MET; 
    if (METnomu > MET) tMET=METnomu; 
    if(tMET < 90.0)  
    { 
      eff=0.0; 
    }  
    else if (tMET < 250.0) 
    { 
      eff = (tMET-50.0)/200.0; 
    } 
    else 
    { 
      eff=1.0; 
    } 
     
    if(rd(engine) > eff) passMET2015=false; 
  } 
 
  if(passMET) 
  { 
    // Use METnomu? 
 
    double tMET=MET; 
    if (METnomu > MET) tMET=METnomu; 
     
    if(tMET < 120.0)  
    { 
      eff=0.0; 
    }  
    else if (tMET < 280.0) 
    { 
      eff = (tMET-120.0)/160.0; 
    } 
    else 
    { 
      eff=1.0; 
    } 
 
    if(rd(engine) > eff) passMET=false; 
  } 

  if(!(Manager()->ApplyCut(passMET,"MET"))) return true; 
  if(!(Manager()->ApplyCut(passMET2015,"MET2015"))) return true; 
//  if(!(Manager()->ApplyCut(passMET2015,"MET2015."))) return true; 
   
  if((!passMET) && (!passMET2015)) {for(auto track: tracks) {delete track;}; return true;} ; 

   
   
   bool onegoodjet = false; 
  bool JetMetDeltaphi=true; 
  bool passJetAngleCuts = true;

  if(event.rec()->jets().size() ==0 ) {for(auto track: tracks) {delete track;}; return true;} ;
 
  std::vector<const RecJetFormat*> Jets;
  for (int i=0; i < event.rec()->jets().size(); i++)
    {
      const RecJetFormat *CurrentJet = &(event.rec()->jets()[i]);
      Jets.push_back(CurrentJet);
    }

  SORTER->sort(Jets);
  std::vector<const RecJetFormat*> GoodJets; 

  if(fabs(Jets[0]->momentum().DeltaPhi(pTmiss)) < 0.5) JetMetDeltaphi = false; 

   filterPhaseSpace(Jets,30.0,2.4);

   for(int m = 0; m< Jets.size() ; m++) 
  {
    
    const RecJetFormat* tjet=Jets[m];
  
  
    if((!onegoodjet) && (fabs(tjet->eta())< 2.4) && (tjet->pt() > 110.0))  
    { 
      onegoodjet=true; 
    } 
    if(tjet->pt() > 30.0)  
    {  
      GoodJets.push_back(tjet); 
 
    } 
    for(int n=m+1; n < Jets.size(); n++) 
    { 
  
      if(fabs(tjet->momentum().DeltaPhi(Jets[n]->momentum())) > 2.5) 
      { 
        passJetAngleCuts = false; 
        break; 
      } 
 
    } 
 
  } 


   if(!(Manager()->ApplyCut(onegoodjet,"OneGoodJet"))) return true; 
  if(!onegoodjet)  {for(auto track: tracks) {delete track;}; return true;} ; 
   
 
     if(!(Manager()->ApplyCut(passJetAngleCuts,"JetAngleCuts"))) return true; 
   if(!(Manager()->ApplyCut(JetMetDeltaphi,"|Delta phi(leading jet, pTmiss)| > 0.5"))) return true; 
  if((!passJetAngleCuts) || (!JetMetDeltaphi))  {for(auto track: tracks) {delete track;}; return true;} ;



  if(tracks.size() < 1)  return true; 
  bool good2015_3=false; 
  bool good2016_3=false; 
  bool good2017_1=false; 
  bool good2017_2=false; 
  bool good2017_3=false; 
  bool good2018A_1=false; 
  bool good2018A_2=false; 
  bool good2018A_3=false; 
  bool good2018B_1=false; 
  bool good2018B_2=false; 
  bool good2018B_3=false; 
 
  bool good2015_3dR=false; 
  bool good2016_3dR=false; 
  bool good2017_1dR=false; 
  bool good2017_2dR=false; 
  bool good2017_3dR=false; 
  bool good2018A_1dR=false; 
  bool good2018A_2dR=false; 
  bool good2018A_3dR=false; 
  bool good2018B_1dR=false; 
  bool good2018B_2dR=false; 
  bool good2018B_3dR=false; 
 
  bool tracketa=false; 
  bool trackpT=false; 
  bool trackfiducial=false; 
  bool trackpixel=false; 
  bool trackinnermiddle=false;  
  bool trackisolation=false; 
  bool trackd0=false; 
  bool trackdz=false; 
  bool trackDRjet=false; 
   
 
  std::vector<charged_track*> tracks2; 

  for(auto chargino: tracks) 
  { 
    if((chargino->abseta() > 2.1) ) {delete chargino; continue;}; 
 
    tracketa=true; 


    if(chargino->Pt() < 55.0) {delete chargino; continue;}; 
    trackpT=true; 
    

    
    double cheta=chargino->abseta(); 

    
    // apply conditions on the track to avoid low efficiency parts of muon chamber/ECAL 
    // muon chamber 
    //if((decayrad > 4020.0) && (decayrad < 7380.0)) 
    //{ 
      if((cheta >0.15) && (cheta < 0.35) ) {delete chargino; continue;}; 
     if((cheta >1.55) && (cheta < 1.85) ) {delete chargino; continue;}; 
    //} 
    // ECAL 
    //if((decayrad > 1290.0) && (decayrad < 1790.0)) 
    //{ 
      if((cheta >1.42) && (cheta < 1.65) ) {delete chargino; continue;}; 
 
    trackfiducial=true; 
 
												
//----------------------------MOD. ISOLATION CALC.----------------------------------------------
    MALorentzVector cprod = chargino->p->ProductionVertex();
    MAfloat32 pv_ssum=0.0, debug_psum=0.0;
    //trackfile << "Event ID: " << eventnumber-1 << endl;

    //Trying to reconstruct results using MCParticles instead of RecTracks, to try to match the isocones number
    for (auto part : event.mc()->particles()){
    	MAbool is_charged = PDG->IsCharged(part.pdgid());
	if (chargino->p->dr(part) < 0.3 && part.statuscode()==1 && !(PHYSICS->Id->IsInvisible(part)) ){ //Maybe should relax statuscode? But would interfere with prop. partons, better leave as is.
		MALorentzVector vprod = part.mothers()[0]->decay_vertex();
		if (vprod.X()<0.01 && vprod.Y()<0.01 && vprod.Z()<0.01) {pv_ssum+=part.momentum().Pt();} //Obs: No definition of displaced in this search, using 0.01mm, according to vertex resolution
		if (part.d0() < 0.2 && part.dz() < 5) {debug_psum+=part.momentum().Pt();}
	}												 //for each coordinate from reconstruction (https://arxiv.org/pdf/1405.6569)
    }
    double isopt = abs(pv_ssum); //Modified from the original to only include contributions coming from the Primary Vertex, no displaced tracks/particles.
    if( (isopt/chargino->Pt()) > 0.05) {delete chargino; continue;}; 
 
    trackisolation=true; 

//--------------------------END MOD. ISOLATION CALC.--------------------------------------------
    												
 
    // Now we have >=1 track with > 3 or 4 pixel hits, actually this means we have to have all 4 closest hits 
 
    if(chargino->_hits.size() < 4) {delete chargino; continue;}; 
 
    bool charginotrackpixel=true; 
     
    for(int i=0; i<4; i++) 
    {  
      if(!chargino->_hits[i]) 
      { 
        charginotrackpixel=false; 
        break; 
      } 
    } 
    if(!charginotrackpixel) {delete chargino; continue;}; 
    trackpixel=true; 
    bool charginoinnermiddle=true; 
    // now check for missing inner/middle hits 
    if(chargino->_hits.size() > 4) 
    { 
      bool foundoutside=false; 
       for(int i=chargino->_hits.size() -1; i>3; i--) 
       { 
         if(foundoutside) 
         { 
           if(!chargino->_hits[i]) 
           { 
             charginoinnermiddle=false; 
             break; 
           } 
 
         } 
         else 
         { 
           if(chargino->_hits[i]) 
           { 
             foundoutside=true; 
           } 
         } 
       }  
 
    } 
 
    if(!charginoinnermiddle) {delete chargino; continue;}; 
 
    trackinnermiddle=true; 

												
//-----------------------------DR CHARGED DAU CUT--------------------------------------------

    //MAbool chdau_iso = true;
    MALorentzVector cdec = chargino->p->mc()->decay_vertex();
    for (auto part : event.mc()->particles()){
	if (abs(part.pdgid())==5){
		MALorentzVector vprod = part.mothers()[0]->decay_vertex(), vdec = part.decay_vertex();
		if (LV_equal(vprod,cdec) && part.pt()>1.0){
			if (chargino->p->dr(part) < 0.2){chargino->chdau_iso = false; break;}
		}
	}
    }

//-----------------------------END DR CHARGED DAU CUT----------------------------------------
												

     
    double absd0=fabs(chargino->p->d0());
    double absdz=fabs(chargino->p->dz());
    
    if(absd0 > 0.2) {delete chargino; continue;}; 
    trackd0=true; 
    if(absdz > 5.0) {delete chargino; continue;}; 
    trackdz=true; 

 
    tracks2.push_back(chargino); 
  } 
 
   if(!(Manager()->ApplyCut(tracketa,">=1 track with |eta| < 2.1"))) return true; 
   if(!(Manager()->ApplyCut(trackpT,">=1 track with p_T > 55"))) return true; 
   if(!(Manager()->ApplyCut(trackfiducial,">=1 track passing fiducial selections"))) return true; 
   if(!(Manager()->ApplyCut(trackisolation,">=1 track with relative isolation"))) return true; 
   if(!(Manager()->ApplyCut(trackpixel,">=1 track with > 3 or 4 pixel hits"))) return true; 
   if(!(Manager()->ApplyCut(trackinnermiddle,">=1 track with no missing inner/middle hits"))) return true; 
   if(!(Manager()->ApplyCut(trackd0,">=1 track with d0 < 0.2 mm"))) return true; 
   if(!(Manager()->ApplyCut(trackdz,">=1 track with dz < 5.0 mm"))) return true; 
 
   
   if(tracks2.size() < 1) {for(auto track: tracks2) {delete track;}; return true;} ; 
 
   if(!(Manager()->ApplyCut(chdau_switch(tracks2),">=1 track with DR(track,daughter)>0.2"))) return true;
 
tracks=FullRemoval(tracks2,GoodJets,0.5); 
 if(!(Manager()->ApplyCut((tracks.size() > 0),">=1 track with DR(track, jet) > 0.5"))) return true; 
 if(!(Manager()->ApplyCut(chdau_switch(tracks),">=1 track with DR(track, jet) > 0.5."))) return true; 
 
if(tracks.size() ==0) return true; 

 tracks=FullRemoval(tracks,event.rec()->electrons(),0.15); 
 if(!(Manager()->ApplyCut((tracks.size() > 0),">=1 track with DR(track, electron) > 0.15"))) return true; 
 if(!(Manager()->ApplyCut(chdau_switch(tracks),">=1 track with DR(track, electron) > 0.15."))) return true; 

if(tracks.size() ==0) return true; 
 
// This cut should remove long-lived charginos which get reconstructed as a muon 
tracks=FullRemoval(tracks,event.rec()->muons(),0.15); 

std::vector<charged_track*> tracks3; 
std::vector<charged_track*> fake_muons; 
  for(auto chargino: tracks) 
  { 
    MALorentzVector vdec = chargino->p->mc()->decay_vertex(); 
    // if it decays outside muon chamber and is not hadronic (i.e. caught in the calo)
    if((!(PHYSICS->Id->IsHadronic(chargino->p->mc()))) &&( (vdec.Pt() > 7000.0) || (fabs(vdec.Pz()) > 11000.0 ))) 
    { 
      fake_muons.push_back(chargino); 
    } 
    else 
    { 
      tracks3.push_back(chargino); 
    } 
  } 
  tracks=FullRemoval(tracks3,fake_muons,0.15); 
  for(auto fmu : fake_muons) 
  { 
    delete fmu; 
  } 
 
 if(!(Manager()->ApplyCut((tracks.size() > 0),">=1 track with DR(track, muon) > 0.15"))) return true; 
 if(!(Manager()->ApplyCut(chdau_switch(tracks),">=1 track with DR(track, muon) > 0.15."))) return true; 
if(tracks.size() ==0) return true; 

 tracks=FullRemoval(tracks,event.rec()->taus(),0.15); 
 if(!(Manager()->ApplyCut((tracks.size() > 0),">=1 track with DR(track, tau) > 0.15"))) return true; 
 if(!(Manager()->ApplyCut(chdau_switch(tracks),">=1 track with DR(track, tau) > 0.15."))) return true; 
 
if(tracks.size() ==0) return true; 
 
bool trackEcalo=false; 
bool trackmissingouter=false; 
 
bool trackEcalo_dR=false;
bool trackmissingouter_dR=false;
 
  for(auto chargino: tracks) 
  { 
   
    bool chtrackEcalo=true;

    /* 
       Now we do the Ecalo, unlike in hackanalysis. The Ecalo is supposed to be total calorimeter
       energy, rather than sumET, but MA5 v1.9 doesn't have that option, so we do what we can.
       I will also include a cut on whether the track is a hadron or electron/tau, under the assumption that 
       in those cases there is a large energy deposit (we should already have thrown it out
       under those cases but ok ...
       
     */
    
    if (PHYSICS->Id->IsHadronic(chargino->p->mc()))
      {
	chtrackEcalo=false;
	continue;
      };

    
    if(chargino->abspid() < 25)
      {
	chtrackEcalo=false;
	continue;
      }
    double isoET=fabs(chargino->p->isolCones()[1].sumET());
    
    if(isoET > 10.0)
      {
	chtrackEcalo=false;
	continue;
      }
    if(!chtrackEcalo) continue;
    
    trackEcalo=true;
    if (chtrackEcalo && chargino->chdau_iso) trackEcalo_dR=true;
    
    double chphi=chargino->phi(); 
    double cheta=chargino->eta(); 
    double acheta=fabs(cheta); 
 
    bool chtrackmissingouter=true; 
    int nhits=(chargino->_hits.size())-1; 
    for(int i=0; i<3; i++)  
       { 
         if(chargino->_hits[nhits-i])  // 3 outermost hits must be missing 
         { 
           chtrackmissingouter=false; 
           break; 
         } 
 
       } 
    if(!chtrackmissingouter) continue; 
    trackmissingouter=true; 
    if(chtrackmissingouter && chargino->chdau_iso) trackmissingouter_dR=true;
    
    bool good2017 = true; 
    bool good2018 = true; 
     
     
 
 
    if((cheta >0.0) && (cheta < 1.42) && (chphi < 3.142) && (chphi > 2.7) ) good2017=false; 
    if((cheta >0.0) && (cheta < 1.42) && (chphi < 0.8) && (chphi > 0.4) ) good2018=false; 
 
    // Now we've already selected tracks that have no missing hits from start to end, so just count layers 
 
    int nlayers=0; 
 
    for(auto hit : chargino->_hits) 
    { 
      if(hit) nlayers++; 
    } 
 
    if(nlayers == 4) 
    { 
      if(good2017) good2017_1=true; 
      if(good2017 && chargino->chdau_iso) good2017_1dR=true;
      
      if(good2018) 
        { 
          good2018A_1=true; 
          good2018B_1=true; 
        } 
      
      if(good2018 && chargino->chdau_iso) 
        { 
          good2018A_1dR=true; 
          good2018B_1dR=true; 
        } 
    } 
    
    else if(nlayers == 5) 
    { 
      if(good2017) good2017_2=true; 
      if(good2017 && chargino->chdau_iso) good2017_2dR=true; 

        if(good2018) 
        { 
          good2018A_2=true; 
          good2018B_2=true; 
        } 
 
        if(good2018 && chargino->chdau_iso) 
        { 
          good2018A_2dR=true; 
          good2018B_2dR=true; 
        } 
 
    } 
    
    else // nlayers > 5 
    { 
        if(nlayers > 7) 
        { 
          good2015_3=true; 
          good2016_3=true; 
        } 
        if(nlayers > 7 && chargino->chdau_iso) 
        { 
          good2015_3dR=true; 
          good2016_3dR=true; 
        } 

        if(good2017) good2017_3=true; 
        
        if(good2017 && chargino->chdau_iso) good2017_3dR=true;
        
        if(good2018) 
        { 
          good2018A_3=true; 
          good2018B_3=true; 
        } 

        if(good2018 && chargino->chdau_iso) 
        { 
          good2018A_3dR=true; 
          good2018B_3dR=true; 
        } 
    } 
 
  } 
 
 if(!(Manager()->ApplyCut(trackEcalo,">=1 track with Ecalo < 10"))) return true; 
 if(!(Manager()->ApplyCut(trackEcalo_dR,">=1 track with Ecalo < 10."))) return true; 
 if(!(Manager()->ApplyCut(trackmissingouter,">=1 track with >=3 missing outer hits"))) return true; 
 if(!(Manager()->ApplyCut(trackmissingouter_dR,">=1 track with >=3 missing outer hits."))) return true; 
 
 
double phipTmiss = pTmiss.Phi();  
bool GoodpTPhi=true; 
if((phipTmiss > -1.6) && (phipTmiss <-0.6)) GoodpTPhi=false; 
 
 if(!(Manager()->ApplyCut(GoodpTPhi,"HEMVeto"))) return true; 
 if(!(Manager()->ApplyCut(GoodpTPhi,"HEMVeto."))) return true; 
 
 
 
 if(!(Manager()->ApplyCut(good2017_1,"4 Layers 2017"))) return true; 
 if(!(Manager()->ApplyCut(good2017_1dR,"4 Layers 2017."))) return true; 
 if(!(Manager()->ApplyCut(good2018A_1,"4 Layers 2018A"))) return true; 
 if(!(Manager()->ApplyCut(good2018A_1dR,"4 Layers 2018A."))) return true; 
 if(!(Manager()->ApplyCut(good2018B_1,"4 Layers 2018B"))) return true; 
 if(!(Manager()->ApplyCut(good2018B_1dR,"4 Layers 2018B."))) return true; 
 
 if(!(Manager()->ApplyCut(good2017_2,"5 Layers 2017"))) return true; 
 if(!(Manager()->ApplyCut(good2017_2dR,"5 Layers 2017."))) return true; 
 if(!(Manager()->ApplyCut(good2018A_2,"5 Layers 2018A"))) return true; 
 if(!(Manager()->ApplyCut(good2018A_2dR,"5 Layers 2018A."))) return true; 
 if(!(Manager()->ApplyCut(good2018B_2,"5 Layers 2018B"))) return true; 
 if(!(Manager()->ApplyCut(good2018B_2dR,"5 Layers 2018B."))) return true; 
 
 if(!(Manager()->ApplyCut(good2015_3,">5 Layers 2015"))) return true; 
 if(!(Manager()->ApplyCut(good2015_3dR,">5 Layers 2015."))) return true; 
 if(!(Manager()->ApplyCut(good2016_3,">5 Layers 2016"))) return true; 
 if(!(Manager()->ApplyCut(good2016_3dR,">5 Layers 2016."))) return true; 
 if(!(Manager()->ApplyCut(good2017_3,">5 Layers 2017"))) return true; 
 if(!(Manager()->ApplyCut(good2017_3dR,">5 Layers 2017."))) return true; 
 if(!(Manager()->ApplyCut(good2018A_3,">5 Layers 2018A"))) return true; 
 if(!(Manager()->ApplyCut(good2018A_3dR,">5 Layers 2018A."))) return true; 
 if(!(Manager()->ApplyCut(good2018B_3,">5 Layers 2018B"))) return true; 
 if(!(Manager()->ApplyCut(good2018B_3dR,">5 Layers 2018B."))) return true; 

for(auto track: tracks) {delete track;}; 


// Now apply reweightings

//Manager()->SetCurrentEventWeight(eventWeight);

 const std::string rwghtname="Reweight ";
 for(auto period : alldataperiods)
   {
      Manager()->SetCurrentEventWeight(myWeight*lumiratios[period]);
      if(!(Manager()->ApplyCut(true,rwghtname+period))) return true;
      if(!(Manager()->ApplyCut(true,rwghtname+period+"."))) return true;
   }

 
 return true; 

}

