#define NUCDEEX_IN_GEANT4_MODE 1
#include "globals.hh"

#ifndef G4NucDeExInterfaceMessenger_hh
#define G4NucDeExInterfaceMessenger_hh

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4String.hh"

class G4NucDeExInterfaceMessenger : public G4UImessenger
{

  public:
    G4NucDeExInterfaceMessenger ();
    ~G4NucDeExInterfaceMessenger ();
    void SetNewValue (G4UIcommand *command, G4String newValues);

  private:
    static const G4String theUIDirectory;
    G4UIdirectory *theNucDeExDirectory;
    G4UIcmdWithAnInteger *verboseCmd;
    G4UIcmdWithAnInteger *seedCmd;
};

#endif

