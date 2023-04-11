#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "constants.hpp"
#include "Model.hpp"
#include "MsRgbModels/GenericMsModel.hpp"
#include "WdCoolingModels/AlthausWdModel.hpp"
#include "WdCoolingModels/MontgomeryWdModel.hpp"
#include "WdCoolingModels/NewMontgomeryWdModel.hpp"
#include "WdCoolingModels/RenedoWdModel.hpp"
#include "WdCoolingModels/WoodWdModel.hpp"
#include "WdAtmosphereModels/BergeronAtmosphereModel.hpp"
#include "WdAtmosphereModels/Bergeron2019AtmosphereModel.hpp"
#include "WdAtmosphereModels/Bergeron2020AtmosphereModel.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::shared_ptr;
using std::string;
using std::vector;

namespace internal
{
    shared_ptr<MsRgbModel> createMsRgbModel(string modelPath)
    {
        return shared_ptr<MsRgbModel>(new GenericMsModel(modelPath));
    }

    shared_ptr<WdCoolingModel> createWdCoolingModel(WdModel model)
    {
        switch (model)
        {
            case WdModel::WOOD:
                return shared_ptr<WoodWdModel>(new WoodWdModel);
            case WdModel::MONTGOMERY:
                return shared_ptr<MontgomeryWdModel>(new MontgomeryWdModel);
            case WdModel::ALTHAUS:
                return shared_ptr<AlthausWdModel>(new AlthausWdModel);
            case WdModel::RENEDO:
                return shared_ptr<RenedoWdModel>(new RenedoWdModel);
            case WdModel::NEW_MONTGOMERY:
                return shared_ptr<NewMontgomeryWdModel>(new NewMontgomeryWdModel);
            default:
                cerr << "***Error: No model found for white dwarf filter set " << static_cast<int>(model) << ".***" << endl;
                cerr << "[Exiting...]" << endl;
                exit (1);
        }
    }

    shared_ptr<WdAtmosphereModel> createWdAtmosphereModel(WdAtmosphereModelSet model)
    {
        switch (model)
        {
            case WdAtmosphereModelSet::BERGERON:
                return shared_ptr<BergeronAtmosphereModel>(new BergeronAtmosphereModel);
            case WdAtmosphereModelSet::BERGERON_2019:
                return shared_ptr<Bergeron2019AtmosphereModel>(new Bergeron2019AtmosphereModel);
            case WdAtmosphereModelSet::BERGERON_2020:
                return shared_ptr<Bergeron2020AtmosphereModel>(new Bergeron2020AtmosphereModel);
            default:
                cerr << "***Error: No model found for white dwarf atmosphere set " << static_cast<int>(model) << ".***" << endl;
                cerr << "[Exiting...]" << endl;
                exit (1);
        }
    }
}

const Model makeModel(const Settings &s)
{
    if (s.verbose)
        cout << "Reading models..." << std::flush;

    Model model( internal::createMsRgbModel(s.mainSequence.modelFile)
               , internal::createWdCoolingModel(s.whiteDwarf.wdModel)
               , internal::createWdAtmosphereModel(s.whiteDwarf.wdAtmosphereModel)
               , s.verbose);

// !!! FIX ME !!!

    model.IFMR = s.whiteDwarf.ifmr;

// END FIX ME

    model.mainSequenceEvol->loadModel(s.files.models);
    model.WDcooling->loadModel(s.files.models);
    model.WDAtmosphere->loadModel(s.files.models);

    if (s.verbose)
        cout << " Done.\n" << endl;

    return model;
}
