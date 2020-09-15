// Copyright (c) 2018, ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: Johannes L. Schoenberger (jsch-at-demuc-dot-de)

#ifndef COLMAP_SRC_OPTIM_LORANSAC_H_
#define COLMAP_SRC_OPTIM_LORANSAC_H_

#include <cfloat>
#include <random>
#include <stdexcept>
#include <vector>

#include "optim/random_sampler.h"
#include "optim/ransac.h"
#include "optim/support_measurement.h"
#include "util/alignment.h"
#include "util/logging.h"

namespace colmap {

// Implementation of LO-RANSAC (Locally Optimized RANSAC).
//
// "Locally Optimized RANSAC" Ondrej Chum, Jiri Matas, Josef Kittler, DAGM 2003.
template <typename Estimator, typename LocalEstimator,
          typename SupportMeasurer = InlierSupportMeasurer,
          typename Sampler = RandomSampler>
class LORANSAC : public RANSAC<Estimator, SupportMeasurer, Sampler> {
 public:
  using typename RANSAC<Estimator, SupportMeasurer, Sampler>::Report;

  explicit LORANSAC(const RANSACOptions& options);

  // Robustly estimate model with RANSAC (RANdom SAmple Consensus).
  //
  // @param X              Independent variables.
  // @param Y              Dependent variables.
  //
  // @return               The report with the results of the estimation.
  Report Estimate(const std::vector<typename Estimator::X_t>& X,
                  const std::vector<typename Estimator::Y_t>& Y);

  // Objects used in RANSAC procedure.
  using RANSAC<Estimator, SupportMeasurer, Sampler>::estimator;
  LocalEstimator local_estimator;
  using RANSAC<Estimator, SupportMeasurer, Sampler>::sampler;
  using RANSAC<Estimator, SupportMeasurer, Sampler>::support_measurer;

 private:
  using RANSAC<Estimator, SupportMeasurer, Sampler>::options_;
};

////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////

template <typename Estimator, typename LocalEstimator, typename SupportMeasurer,
          typename Sampler>
LORANSAC<Estimator, LocalEstimator, SupportMeasurer, Sampler>::LORANSAC(
    const RANSACOptions& options)
    : RANSAC<Estimator, SupportMeasurer, Sampler>(options) {}

template <typename Estimator, typename LocalEstimator, typename SupportMeasurer,
          typename Sampler>
typename LORANSAC<Estimator, LocalEstimator, SupportMeasurer, Sampler>::Report
LORANSAC<Estimator, LocalEstimator, SupportMeasurer, Sampler>::Estimate(
    const std::vector<typename Estimator::X_t>& X,
    const std::vector<typename Estimator::Y_t>& Y) {
  CHECK_EQ(X.size(), Y.size());

  const size_t num_samples = X.size();

  typename RANSAC<Estimator, SupportMeasurer, Sampler>::Report report;
  report.success = false;
  report.num_trials = 0;

  if (num_samples < Estimator::kMinNumSamples) {
    return report;
  }

  typename SupportMeasurer::Support best_support;
  typename Estimator::M_t best_model;
  bool best_model_is_local = false;

  bool abort = false;

  const double max_residual = options_.max_error * options_.max_error;

  std::vector<double> residuals(num_samples);

  std::vector<typename LocalEstimator::X_t> X_inlier;
  std::vector<typename LocalEstimator::Y_t> Y_inlier;

  std::vector<typename Estimator::X_t> X_rand(Estimator::kMinNumSamples);
  std::vector<typename Estimator::Y_t> Y_rand(Estimator::kMinNumSamples);

  sampler.Initialize(num_samples);

  size_t max_num_trials = options_.max_num_trials;
  max_num_trials = std::min<size_t>(max_num_trials, sampler.MaxNumSamples());
  size_t dyn_max_num_trials = max_num_trials;

  for (report.num_trials = 0; report.num_trials < max_num_trials;
       ++report.num_trials) {
    //std::cout<<"////////report.num_trials: "<<report.num_trials<<"\n";
    if (abort) {
      report.num_trials += 1;
      break;
    }

    sampler.SampleXY(X, Y, &X_rand, &Y_rand);

    // Estimate model for current subset.
    const std::vector<typename Estimator::M_t> sample_models =
        estimator.Estimate(X_rand, Y_rand);

    // Iterate through all estimated models
    for (const auto& sample_model : sample_models) {
      estimator.Residuals(X, Y, sample_model, &residuals);
      CHECK_EQ(residuals.size(), X.size());

      const auto support = support_measurer.Evaluate(residuals, max_residual);
      //std::cout<<"////////support.num_inliers: "<<support.num_inliers<<"\n";
      
      // Do local optimization if better than all previous subsets.
      if (support_measurer.Compare(support, best_support)) {
        best_support = support;
        best_model = sample_model;
        best_model_is_local = false;
        
        // Estimate locally optimized model from inliers.
        if (support.num_inliers > Estimator::kMinNumSamples &&
            support.num_inliers >= LocalEstimator::kMinNumSamples) {
        //std::cout<<"////////localransac: \n";
          X_inlier.clear();
          Y_inlier.clear();
          X_inlier.reserve(support.num_inliers);
          Y_inlier.reserve(support.num_inliers);
          for (size_t i = 0; i < residuals.size(); ++i) {
            if (residuals[i] <= max_residual) {
              X_inlier.push_back(X[i]);
              Y_inlier.push_back(Y[i]);
            }
          }

          const std::vector<typename LocalEstimator::M_t> local_models =
              local_estimator.Estimate(X_inlier, Y_inlier);

          for (const auto& local_model : local_models) {
            local_estimator.Residuals(X, Y, local_model, &residuals);
            CHECK_EQ(residuals.size(), X.size());

            const auto local_support =
                support_measurer.Evaluate(residuals, max_residual);
            //std::cout<<"////////local_support.num_inliers: "<<local_support.num_inliers<<"\n";
            // Check if non-locally optimized model is better.
            if (support_measurer.Compare(local_support, best_support)) {
              //std::cout<<"////////local_support.is better\n";
              best_support = local_support;
              best_model = local_model;
              best_model_is_local = true;
            }
          }
        }
        //std::cout<<"/////dyn_max_num_trials: "<<dyn_max_num_trials<<"\n";

        dyn_max_num_trials =
            RANSAC<Estimator, SupportMeasurer, Sampler>::ComputeNumTrials(
                best_support.num_inliers, num_samples, options_.confidence,
                options_.dyn_num_trials_multiplier);
        //std::cout<<"/////recalculated dyn_max_num_trials: "<<dyn_max_num_trials<<"\n";
      }

      if (report.num_trials >= dyn_max_num_trials &&
          report.num_trials >= options_.min_num_trials) {
        // std::cout<<"/////abort:  report.num_trials("<<report.num_trials<<") >="
        //          <<"dyn_max_num_trials("<<dyn_max_num_trials<<")\n";
                 
        abort = true;
        break;
      }
    }
  }

  report.support = best_support;
  report.model = best_model;

  // No valid model was found
  if (report.support.num_inliers < estimator.kMinNumSamples) {
    return report;
  }

  report.success = true;

  // Determine inlier mask. Note that this calculates the residuals for the
  // best model twice, but saves to copy and fill the inlier mask for each
  // evaluated model. Some benchmarking revealed that this approach is faster.

  if (best_model_is_local) {
    local_estimator.Residuals(X, Y, report.model, &residuals);
  } else {
    estimator.Residuals(X, Y, report.model, &residuals);
  }

  CHECK_EQ(residuals.size(), X.size());

  report.inlier_mask.resize(num_samples);
  for (size_t i = 0; i < residuals.size(); ++i) {
    if (residuals[i] <= max_residual) {
      report.inlier_mask[i] = true;
    } else {
      report.inlier_mask[i] = false;
    }
  }

  return report;
}

}  // namespace colmap

#endif  // COLMAP_SRC_OPTIM_LORANSAC_H_
