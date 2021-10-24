//
//  Inner proton absorber.  Currently a thin cylinder
//
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
#include "TrackToy/Detector/CylindricalShell.hh"
namespace TrackToy {
  class IPA {
    public:
      enum IPAType { unknown=-1, cylinder=1, propeller };
      IPA(MatEnv::MatDBInfo const& matdbinfo,std::string const& file);
      auto const& cyl() const { return cyl_; }
      auto const& material() const { return *mat_; }
    private:
      IPAType type_;
      CylindricalShell cyl_;
      const MatEnv::DetMaterial* mat_;

  };
}


