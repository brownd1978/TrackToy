//
//  Inner proton absorber.  Currently a thin cylinder
//
#include "TrackToy/Detector/CylindricalShell.hh"
namespace TrackToy {
  class IPA {
    public:
      enum IPAType { unknown=-1, cylinder=1, propeller };
      IPA(std::string const& file);
      auto const& cylinder() const { return cyl_; }
    private:
      IPAType type_;
      CylindricalShell cyl_;

  };
}


