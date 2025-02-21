#pragma once

#include "nusystematics/utility/exceptions.hh"
#include "systematicstools/interface/ISystProviderTool.hh"

#include "fhiclcpp/ParameterSet.h"

// GENIE
#include "Framework/EventGen/EventRecord.h"
// Extra includes needed for CheckTune()
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/XSecSplineList.h"

namespace nusyst {

class IGENIESystProvider_tool : public systtools::ISystProviderTool {
protected:
  // Based on GENIEHelper::FindTune()
  // TODO: reduce code duplication here
  // -- S. Gardiner, 20 December 2018
  void CheckTune(const std::string &tune_name) {

    std::string fhicl_tune_name = tune_name;

    // The default tune name is ${GENIE_XSEC_TUNE}, which
    // should be converted into the value of the corresponding
    // enviornment variable, as is done below.
    if (fhicl_tune_name.front() == '$') {
      // need to remove ${}'s
      std::string tuneEnvVar = fhicl_tune_name;
      std::string rmchars("$(){} ");
      // std::remove_if removes characters in [first,last) that are found
      //   within the rmchars string. It returns returns a past-the-end
      //   iterator for the new end of the range [funky!]
      // std::string::erase actually trims the string
      tuneEnvVar.erase(std::remove_if(tuneEnvVar.begin(), tuneEnvVar.end(),
                                      [&rmchars](const char &c) -> bool {
                                        return rmchars.find(c) !=
                                               std::string::npos;
                                      }),
                       tuneEnvVar.end());

      const char *tune = std::getenv(tuneEnvVar.c_str());
      if (tune) {
        fhicl_tune_name = std::string(tune);
      } else {
        throw systtools::invalid_ToolConfigurationFHiCL()
            << "can't resolve TuneName: " << fhicl_tune_name;
      }
    }

    // If the XSecSplineList returns a non-empty string as the current tune
    // name, then genie::RunOpt::BuildTune() has already been called.
    std::string current_tune = genie::XSecSplineList::Instance()->CurrentTune();
    if (current_tune.empty()) {
      // Constructor automatically calls grunopt->Init();
      genie::RunOpt *grunopt = genie::RunOpt::Instance();
      grunopt->SetTuneName(fhicl_tune_name);
      grunopt->BuildTune();
    } else {
      // It has already been built, so just check consistency
      if (fhicl_tune_name != current_tune) {
        throw systtools::invalid_ToolConfigurationFHiCL()
            << "Requested GENIE tune \"" << fhicl_tune_name
            << "\" does not match previously built tune \"" << current_tune
            << '\"';
      }
    }
  }

public:
  IGENIESystProvider_tool(fhicl::ParameterSet const &ps)
      : ISystProviderTool(ps), fGENIEModuleLabel(ps.get<std::string>(
                                   "genie_module_label", "generator")) {

    std::string tune_name =
        ps.get<std::string>("TuneName", "${GENIE_XSEC_TUNE}");
    this->CheckTune(tune_name);
  }

  NEW_SYSTTOOLS_EXCEPT(invalid_response);

  /// Calculates configured response for a given GHep record
  virtual systtools::event_unit_response_t
  GetEventResponse(genie::EventRecord const &) = 0;

  /// Calculates response, allowing parameter values to be overridden locally
  systtools::event_unit_response_t
  GetEventResponse(genie::EventRecord const & evt,
    std::vector<std::pair<systtools::paramId_t, std::vector<double>>> const & paramVals)
  {
    systtools::SystMetaData params = this->GetSystMetaData();
    for (const auto & paramList  : paramVals)
    {

      systtools::paramId_t param = paramList.first;
      auto it_Vars = std::find_if(params.begin(), params.end(),
                                                   [param](const systtools::SystParamHeader & p) { return p.systParamId == param; });

      if (it_Vars == params.end())
        throw systtools::parameter_Id_not_handled() << "Parameter " << std::to_string(param)
                                                    << " is not configured, can't calculate response for it";

      this->OverrideVariations(param, paramList.second);
    }

    return GetEventResponse(evt);
  }

  /// Calculates configured response for a given vector of GHep record
  std::unique_ptr<systtools::EventResponse>
  GetEventResponses(std::vector<std::unique_ptr<genie::EventRecord>> const &gheps){

    std::unique_ptr<systtools::EventResponse> er =
        std::make_unique<systtools::EventResponse>();

    for (size_t eu_it = 0; eu_it < gheps.size(); ++eu_it) {
      er->push_back(GetEventResponse(*gheps[eu_it]));
    }
    return er;

  };

  systtools::event_unit_response_w_cv_t
  GetEventVariationAndCVResponse(genie::EventRecord const &GenieGHep) {
    systtools::event_unit_response_w_cv_t responseandCV;

    systtools::event_unit_response_t prov_response =
        GetEventResponse(GenieGHep);

    // Foreach param
    for (systtools::ParamResponses &pr : prov_response) {
      // Get CV resp
      systtools::SystParamHeader const &hdr =
          GetParam(GetSystMetaData(), pr.pid);

      // If not a correction dial, responses and paramVariations should have same size
      if ( !hdr.isCorrection &&
           (pr.responses.size() != hdr.paramVariations.size())
      ) {
        throw invalid_response()
            << "[ERROR]: Parameter: " << hdr.prettyName << ", with "
            << hdr.paramVariations.size() << " parameter variations, returned "
            << pr.responses.size() << " responses.";
      }
      // make sure correction dial has zero paramVariations.size()
      if( hdr.isCorrection && hdr.paramVariations.size() != 0 ){
        throw invalid_response()
            << "[ERROR]: Parameter: " << hdr.prettyName << " is a correction but has non-zero parameter variations ("
            << hdr.paramVariations.size() << ").";
      }
      // make sure correction dial has exactly one response
      if( hdr.isCorrection && pr.responses.size() != 1 ){
        throw invalid_response()
            << "[ERROR]: Parameter: " << hdr.prettyName << " is a correction and should have single response, but got"
            << pr.responses.size() << " responses.";
      }

      double CVResp = hdr.isWeightSystematicVariation ? 1 : 0;
      size_t NVars = hdr.paramVariations.size(); // note: NVars is zero for correction dial

      // If CV is different from default, find it from paramVariations and get the CV weight,
      // then divide all the weights by this CV weight.
      // Analyzers should apply the CV weight first and then multiply each response
      if (hdr.centralParamValue != systtools::kDefaultDouble) {
        // note: NVars is zero for correction dial
        for (size_t idx = 0; idx < NVars; ++idx) {
          if (fabs(hdr.centralParamValue - hdr.paramVariations[idx]) <=
              std::numeric_limits<float>::epsilon()) {
            CVResp = pr.responses[idx];
            break;
          }
        }
        // if we didn't find it, the CVResp stays as 1/0 depending on whether it
        // is a weight or not.
        for (size_t idx = 0; idx < NVars; ++idx) {
          if (hdr.isWeightSystematicVariation) {
            pr.responses[idx] /= CVResp; // divide the responses by CV weight
          } else {
            pr.responses[idx] -= CVResp;
          }
        }
        // For a correction dial, we have NVars=0, so manually update CVResp and responses
        if( hdr.isCorrection ){
          CVResp = pr.responses[0];
          pr.responses[0] = 1.;
        }
      }

      responseandCV.push_back({pr.pid, CVResp, pr.responses});
    } // end for parameter response

    return responseandCV;
  }

  /// Calculates the response to a single parameter for a given GHep record
  virtual systtools::event_unit_response_t
  GetEventResponse(genie::EventRecord const &, systtools::paramId_t) {
    throw systtools::ISystProviderTool_method_unimplemented()
        << "[ERROR]: " << GetFullyQualifiedName()
        << " does not implement systtools::event_unit_response_t "
           "GetEventResponse(genie::EventRecord &, systtools::paramId_t).";
  }

  /// Calculates the multiplicatively combined responses for a given set of
  /// parameter--value pairs.
  ///
  /// \note This convenience method should only be used for weight responses.
  virtual double GetEventWeightResponse(genie::EventRecord const &,
                                        systtools::param_value_list_t const &) {
    throw systtools::ISystProviderTool_method_unimplemented()
        << "[ERROR]: " << GetFullyQualifiedName()
        << " does not implement double "
           "GetEventWeightResponse(genie::EventRecord "
           "&,systtools::param_value_list_t const &).";
  }

  std::string fGENIEModuleLabel;
};
} // namespace nusyst
