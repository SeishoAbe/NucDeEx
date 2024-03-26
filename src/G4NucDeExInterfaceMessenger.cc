#define NUCDEEX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4NucDeExInterfaceMessenger.hh"
#include "NucDeExUtils.hh"
#include "NucDeExRandom.hh"
#include "G4SystemOfUnits.hh"
#include <sstream>

const G4String G4NucDeExInterfaceMessenger::theUIDirectory = "/process/had/nucdeex/";

G4NucDeExInterfaceMessenger::G4NucDeExInterfaceMessenger(G4NucDeExInterface* anInterface):
  theNucDeExInterface(anInterface)
{
  // Create a directory for the INCL++ commands
  theNucDeExDirectory = new G4UIdirectory(theUIDirectory);
  theNucDeExDirectory->SetGuidance("Parameters for the NucDeEx model");

  verboseCmd = new G4UIcmdWithAnInteger((theUIDirectory + "verbose").data(),this);
  verboseCmd->SetGuidance("Set verbosity of NucDeEx");
  verboseCmd->SetParameterName("NucDeExVerbose",true);
  verboseCmd->SetDefaultValue(0);
  verboseCmd->SetRange("NucDeExVerbose>=0");
  verboseCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  seedCmd = new G4UIcmdWithAnInteger((theUIDirectory + "seed").data(),this);
  seedCmd->SetGuidance("Set seed of NucDeEx");
  seedCmd->SetParameterName("NucDeExseed",true);
  seedCmd->SetDefaultValue(1);
  seedCmd->SetRange("NucDeExseed>=0");
  seedCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

G4NucDeExInterfaceMessenger::~G4NucDeExInterfaceMessenger() {
  delete theNucDeExDirectory;
  delete verboseCmd;
  delete seedCmd;
}

void G4NucDeExInterfaceMessenger::SetNewValue(G4UIcommand *command, G4String newValues) {
  if(command==verboseCmd) {
    const G4int parameter = verboseCmd->GetNewIntValue(newValues);
    NucDeEx::Utils::fVerbose=parameter;
  }else if(command==seedCmd) {
    const G4int parameter = seedCmd->GetNewIntValue(newValues);
    NucDeEx::Random::SetSeed(parameter);
  }
}
