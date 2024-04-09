#ifndef analysis_cms_exo_19_010_h
#define analysis_cms_exo_19_010_h

#include "SampleAnalyzer/Process/Analyzer/AnalyzerBase.h"
								
#include "SampleAnalyzer/Commons/Service/PDGService.h"
								
//#include "new_smearer_reco.h"
//#include "new_tagger.h"
#include <random>
namespace MA5
{
class cms_exo_19_010 : public AnalyzerBase
{
  INIT_ANALYSIS(cms_exo_19_010,"cms_exo_19_010")

 public:
  virtual bool Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters);
  virtual void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
  virtual bool Execute(SampleFormat& sample, const EventFormat& event);

 private:
  std::random_device m_randomdevice;
  std::mt19937 engine;

  std::map<std::string, double> lumiratios;
  const std::vector<std::string> alldataperiods = {"2015", "2016A", "2016B", "2017", "2018A", "2018B"};
};
}

template<typename T1, typename T2> std::vector<T1*>
  FullRemoval(std::vector<T1*> &v1, const std::vector<T2> &v2,
  const double &drmin)
{
  // Determining with objects should be removed, and free memory
  
  std::vector<bool> mask(v1.size(),false);

  for (unsigned int j=0;j<v1.size();j++)
  {
    for (unsigned int i=0; !mask[j] && i<v2.size();i++)
        {
            if (v2[i].momentum().DeltaR(v1[j]->momentum()) < drmin) 
	            {
	                mask[j]=true;
	                //break;
	            }
        }
    };
  // Building the cleaned container
  std::vector<T1*> cleaned_v1;
  for (unsigned int i=0;i<v1.size();i++)
  {
    if (!mask[i]) { 
      cleaned_v1.push_back(v1[i]);
      }
      else
      {
        delete v1[i];
      }

  }

  return cleaned_v1;
};



template<typename T1, typename T2> std::vector<T1*>
  FullRemoval(std::vector<T1*> &v1, std::vector<T2*> &v2,
  const double &drmin)
{
  // Determining with objects should be removed, and free memory
  
  std::vector<bool> mask(v1.size(),false);

  for (unsigned int j=0;j<v1.size();j++)
  {
    for (unsigned int i=0; !mask[j] && i<v2.size();i++)
        {
            if (v2[i]->momentum().DeltaR(v1[j]->momentum()) < drmin) 
	            {
	                mask[j]=true;
	                //break;
	            }
        }
    };
  // Building the cleaned container
  std::vector<T1*> cleaned_v1;
  for (unsigned int i=0;i<v1.size();i++)
  {
    if (!mask[i]) { 
      cleaned_v1.push_back(v1[i]);
      }
      else
      {
        delete v1[i];
      }

  }

  return cleaned_v1;
};




template<typename T1, typename T2> std::vector<T1*>
  Removal(std::vector<T1*> &v1, std::vector<T2*> &v2,
  const double &drmin)
{
  // Determining with objects should be removed
  
  std::vector<bool> mask(v1.size(),false);

  for (unsigned int j=0;j<v1.size();j++)
  {
    for (unsigned int i=0; !mask[j] && i<v2.size();i++)
        {
            if (v2[i]->momentum().deltaR(v1[j]->momentum()) < drmin) 
	            {
	                mask[j]=true;
	                //break;
	            }
        }
    };
  // Building the cleaned container
  std::vector<T1*> cleaned_v1;
  for (unsigned int i=0;i<v1.size();i++)
    if (!mask[i]) cleaned_v1.push_back(v1[i]);

  return cleaned_v1;
};


template<typename T1> std::vector<T1*> SelfRemoval(std::vector<T1*> &v1, const double &drmin)
{
  // Determining with objects should be removed -- keep the higher pT version. This only applies for electrons anyway
  
  std::vector<bool> mask(v1.size(),false);

      double tdr;
      for (unsigned int j=0;j<v1.size();j++){
	      for (unsigned int i=j+1;i<v1.size();i++) {

	        tdr=v1[i]->momentum().DeltaR(v1[j]->momentum());
	        if ((tdr < drmin)) 
		      {
           if(v1[i]->momentum().Pt() < v1[j]->momentum().Pt())
            {
		          mask[i]=true;
             }
            else
            {
              mask[j]=true;
            }
		      }

	      }
      };
   

  // Building the cleaned container
  std::vector<T1*> cleaned_v1;
  for (unsigned int i=0;i<v1.size();i++)
    if (!mask[i]) cleaned_v1.push_back(v1[i]);

  return cleaned_v1;
};

template<typename T1> void filterPhaseSpace(std::vector<T1*> &vec_t1, const double &pTmin, const double &absEtaMax)
{
	// filters the particles/jets in a vector by pT and eta.

  
  auto it = vec_t1.begin();
  while(it != vec_t1.end())
  {
	  //T1* t1 = vec_t1[it];
	  //if( (t1->pT() > pTmin) || (t1->abseta() > absEtaMax))
    if( ((*it)->pt() < pTmin ) || (fabs((*it)->eta()) > absEtaMax))
	  {
       it = vec_t1.erase(it);
        }
    else
      {
	it++;
      }
    
  }
  
}



#endif
