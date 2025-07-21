#include "AngleCalibration.hpp"

#include "aare/logger.hpp"

namespace angcal {

AngleCalibration::AngleCalibration(
    std::shared_ptr<MythenDetectorSpecifications> mythen_detector_,
    std::shared_ptr<FlatField> flat_field_,
    const std::filesystem::path &file_path_,
    std::optional<std::shared_ptr<MythenFileReader>> mythen_file_reader_,
    std::optional<std::shared_ptr<SimpleFileInterface>> custom_file_ptr_)
    : mythen_detector(mythen_detector_), flat_field(flat_field_) {

    if (mythen_file_reader_.has_value()) {
        mythen_file_reader = mythen_file_reader_.value();
    } else {
        mythen_file_reader = std::make_shared<MythenFileReader>(file_path_);
    }

    if (custom_file_ptr_.has_value()) {
        custom_file_ptr = custom_file_ptr_.value();
    }

    // calculate inverse normalized flatfield
    flat_field->inverse_normalized_flatfield();

    DGparameters = DGParameters(mythen_detector->max_modules());

    /*
    num_bins = std::floor(mythen_detector->max_angle() / histogram_bin_width)
    -
    std::floor(mythen_detector->min_angle() / histogram_bin_width) +
        1; // TODO only works if negative
           // and positive angle
    */
    num_bins = base_peak_roi * 2 + 1;
}

void AngleCalibration::set_histogram_bin_width(double bin_width) {
    histogram_bin_width = bin_width;

    /*
    num_bins = std::floor(mythen_detector->max_angle() / histogram_bin_width) -
               std::floor(mythen_detector->min_angle() / histogram_bin_width) +
               1; // TODO only works if negative
                  // and positive angle
    */
}

double AngleCalibration::get_histogram_bin_width() const {
    return histogram_bin_width;
}

ssize_t AngleCalibration::get_new_num_bins() const { return num_bins; }

const DGParameters &AngleCalibration::get_DGparameters() const {
    return DGparameters;
}

NDArray<double, 1> AngleCalibration::get_new_photon_counts() const {
    return new_photon_counts;
}

NDArray<double, 1> AngleCalibration::get_new_statistical_errors() const {
    return new_photon_count_errors;
}

void AngleCalibration::read_initial_calibration_from_file(
    const std::string &filename) {

    custom_file_ptr->open(filename);
    custom_file_ptr->read_into(DGparameters.parameters.buffer(), 8);
}

EEParameters AngleCalibration::convert_to_EE_parameters() const {

    EEParameters EEparameters(DGparameters.parameters.shape(0));

    for (ssize_t i = 0; i < DGparameters.parameters.shape(0); ++i) {
        auto [module_center_distance, normal_distance, angle] =
            convert_to_EE_parameters(i);
        EEparameters.normal_distances(i) = normal_distance;
        EEparameters.module_center_distances(i) = module_center_distance;
        EEparameters.angles(i) = angle;
    }

    return EEparameters;
}

std::tuple<double, double, double>
AngleCalibration::convert_to_EE_parameters(const size_t module_index) const {
    return convert_to_EE_parameters(DGparameters.centers(module_index),
                                    DGparameters.conversions(module_index),
                                    DGparameters.offsets(module_index));
}

std::tuple<double, double, double> AngleCalibration::convert_to_EE_parameters(
    const double center, const double conversion, const double offset) const {
    const double module_center_distance =
        center * MythenDetectorSpecifications::pitch();
    const double normal_distance =
        MythenDetectorSpecifications::pitch() / std::abs(conversion);
    const double angle = offset + 180.0 / M_PI * center * std::abs(conversion);

    return std::make_tuple(module_center_distance, normal_distance, angle);
}

size_t AngleCalibration::global_to_local_strip_index_conversion(
    const size_t global_strip_index) const {
    const size_t module_index =
        global_strip_index / MythenDetectorSpecifications::strips_per_module();
    // local strip index in module
    size_t local_strip_index =
        global_strip_index -
        module_index * MythenDetectorSpecifications::strips_per_module();
    // switch if indexing is in clock-wise direction
    local_strip_index =
        std::signbit(DGparameters.conversions(module_index))
            ? MythenDetectorSpecifications::strips_per_module() - 1 -
                  local_strip_index
            : local_strip_index;

    return local_strip_index;
}

/*
parameters
AngleCalibration::convert_to_BC_parameters() {}
*/

double AngleCalibration::diffraction_angle_from_DG_parameters(
    const double center, const double conversion, const double offset,
    const size_t strip_index, const double distance_to_strip) const {

    return offset + 180.0 / M_PI *
                        (center * std::abs(conversion) -
                         atan((center - (strip_index + distance_to_strip)) *
                              std::abs(conversion)));
}

double AngleCalibration::diffraction_angle_from_EE_parameters(
    const double module_center_distance, const double normal_distance,
    const double angle, const size_t strip_index,
    const double distance_to_strip) const {

    return angle - 180.0 / M_PI *
                       atan((module_center_distance -
                             MythenDetectorSpecifications::pitch() *
                                 (strip_index + distance_to_strip)) /
                            normal_distance); // TODO: why is it minus
                                              // is it defined counter
                                              // clockwise? thought
                                              // should have a flipped
                                              // sign
}

double AngleCalibration::angular_strip_width_from_DG_parameters(
    const double center, const double conversion, const double offset,
    const size_t local_strip_index) const {

    return std::abs(diffraction_angle_from_DG_parameters(
                        center, conversion, offset, local_strip_index, -0.5) -
                    diffraction_angle_from_DG_parameters(
                        center, conversion, offset, local_strip_index, 0.5));
}

double AngleCalibration::angular_strip_width_from_EE_parameters(
    const double module_center_distance, const double normal_distance,
    const double angle, const size_t local_strip_index) const {

    return std::abs(diffraction_angle_from_EE_parameters(
                        module_center_distance, normal_distance, angle,
                        local_strip_index, -0.5) -
                    diffraction_angle_from_EE_parameters(
                        module_center_distance, normal_distance, angle,
                        local_strip_index, 0.5));

    // TODO: again not sure about division order - taking abs anyway
}

double
AngleCalibration::similarity_criterion(const NDView<double, 1> S0,
                                       const NDView<double, 1> S1,
                                       const NDView<double, 1> S2) const {

    double similarity_criterion = 0;
    size_t num_runs = 0;
    for (ssize_t bin_index = 0; bin_index < S0.size(); ++bin_index) {

        num_runs +=
            (S0(bin_index) >
             std::numeric_limits<double>::epsilon()); // doesnt really reflect
                                                      // number of runs and
                                                      // number of bins

        double weighted_average =
            1. / S0(bin_index); // photon variance over each run
        double goodness_of_fit =
            S2(bin_index) -
            S1(bin_index) * S1(bin_index) *
                weighted_average; // calculates chi value for optimal parameter
                                  // a over all runs chi_bin = (a-
                                  // photon_count)*photon_variance
        similarity_criterion +=
            goodness_of_fit *
            weighted_average; // TODO in antonios code only goodness of fit is
                              // added not averaged !!
    }

    return similarity_criterion /
           std::max(num_runs, static_cast<size_t>(
                                  1)); // should actually only divide by runs
}

double AngleCalibration::calculate_similarity_of_peaks(
    const size_t module_index) const {
    // used to calculate similarity criterion between peaks of different
    // acquisition S_index = sum_i^num_runs
    // photon_count^index*photon_variance
    NDArray<double, 1> S2(std::array<ssize_t, 1>{num_bins}, 0.0);
    NDArray<double, 1> S1(std::array<ssize_t, 1>{num_bins}, 0.0);
    NDArray<double, 1> S0(std::array<ssize_t, 1>{num_bins}, 0.0);

    double base_peak_hwid = 0.18; // TODO: not sure what this is - should it be
                                  // configurable

    for (const auto &file : file_list) {
        MythenFrame frame = mythen_file_reader->read_frame(file);

        //++sum_photon_counts; // dont know if intented - just used to skew
        // the distribution?

        // check if base_peak_angle is within the range of angles for this
        // module
        double left_module_boundary_angle =
            diffraction_angle_from_DG_parameters(
                DGparameters.centers(module_index),
                DGparameters.conversions(module_index),
                DGparameters.offsets(module_index), 0);
        double right_module_boundary_angle =
            diffraction_angle_from_DG_parameters(
                DGparameters.centers(module_index),
                DGparameters.conversions(module_index),
                DGparameters.offsets(module_index),
                mythen_detector->strips_per_module());

        // TODO: in antonios code only bloffset and angle is added why? -
        // maybe we can directly add it or is it used later?- careful
        left_module_boundary_angle +=
            (frame.detector_angle + mythen_detector->dtt0() +
             mythen_detector->bloffset()); // Antonio didnt add dtt0

        right_module_boundary_angle +=
            (frame.detector_angle + mythen_detector->dtt0() +
             mythen_detector->bloffset());

        if (base_peak_angle + base_peak_hwid < left_module_boundary_angle ||
            base_peak_angle > right_module_boundary_angle)
            continue; // skip module if base peak angle is not in range

        // we should have histograms per run!!!
        /*
        for (size_t strip_index =
                 module_index * mythen_detector->strips_per_module();
             strip_index <
             (module_index + 1) * mythen_detector->strips_per_module();
             ++strip_index) {

            // TODO: supposed to normalize this thing with
            // total_number_of_photon_count per run and module multiplied
            // with scale factor
            // /(total_number_of_photon_counts_per_run_and_module *
            // scale_factor)
            ++num_runs;
            total_photon_counts(module_index) +=
                frame.photon_counts(strip_index);
        }
        */

        // calculates flatfield normalized photon counts and
        // photon_count_variance for ROI around base_peak and redistributes
        // to fixed angle width bins
        auto [fixed_angle_width_bins_photon_counts,
              fixed_angle_width_bins_photon_count_variance] =
            redistribute_photon_counts_to_fixed_angle_bins(
                module_index, frame, S0.view(), S1.view(), S2.view());
    }

    return similarity_criterion(S0.view(), S1.view(), S2.view());

    /*
    sum_photon_counts += std::accumulate(total_photon_counts.begin(),
                                         total_photon_counts.end(), 0.0);

    // in theory this should happen - but need to store seperately for eacg
    // run

    for (size_t strip_index = 0;
         strip_index < mythen_detector->num_strips(); ++strip_index) {

        size_t module_index =
            strip_index % mythen_detector->strips_per_module();

        flatfield_normalized_photon_counts(strip_index) /=
    (total_photon_counts(module_index) * scale_factor);
        flatfield_normalized_photon_counts_error(strip_index) /=
            (total_photon_counts(module_index) * scale_factor);
    }
    */
}

// TODO: maybe have a function where we calculate the average photon counts -
// multiple loops - or directly store normalized data somewhere instead of
// computing on the fly

void AngleCalibration::calibrate(const std::vector<std::string> &file_list_,
                                 const double base_peak_angle_) {

    file_list = file_list_;

    base_peak_angle = base_peak_angle_;

    for (size_t module_index = 0; module_index < mythen_detector->max_modules();
         ++module_index) {

        LOG(aare::TLogLevel::logINFO)
            << "starting calibration for module " << module_index;

        optimization_algorithm(module_index);
    }
}

// TODO go over function again - in particular peak position and detector
// position
// TODO: maybe store as one 2d array - better cache usage or struct with named
// access functions
std::tuple<NDArray<double, 1>, NDArray<double, 1>>
AngleCalibration::redistribute_photon_counts_to_fixed_angle_bins(
    const size_t module_index, const MythenFrame &frame, NDView<double, 1> S0,
    NDView<double, 1> S1, NDView<double, 1> S2) const {

    NDArray<double, 1> new_fixed_angle_width_bins_photon_variance =
        NDArray<double, 1>(std::array<ssize_t, 1>{num_bins}, 0.0);

    NDArray<double, 1> new_fixed_angle_width_bins_photon_counts =
        NDArray<double, 1>(std::array<ssize_t, 1>{num_bins}, 0.0);

    size_t base_peak_as_bin_index = static_cast<size_t>(
        base_peak_angle /
        histogram_bin_width); // TODO: in antonios code actually rounded to
                              // nearest integer
    double left_boundary_roi_base_peak =
        (base_peak_as_bin_index - base_peak_roi - 0.5) *
        histogram_bin_width; // in degrees //TODO: Antonio subtracts the
                             // detector angle but why?
    double right_boundary_roi_base_peak =
        (base_peak_as_bin_index + base_peak_roi + 0.5) *
        histogram_bin_width; // in degrees

    for (size_t strip_index =
             module_index * mythen_detector->strips_per_module();
         strip_index <
         (module_index + 1) * mythen_detector->strips_per_module();
         ++strip_index) {

        size_t local_strip_index =
            global_to_local_strip_index_conversion(strip_index);
        double left_strip_boundary_angle = diffraction_angle_from_DG_parameters(
            DGparameters.centers(module_index),
            DGparameters.conversions(module_index),
            DGparameters.offsets(module_index), local_strip_index, 0.5);
        double right_strip_boundary_angle =
            diffraction_angle_from_DG_parameters(
                DGparameters.centers(module_index),
                DGparameters.conversions(module_index),
                DGparameters.offsets(module_index), local_strip_index, -0.5);

        if (mythen_detector->get_bad_channels()[strip_index]) {
            continue; // skip bad channels
        }
        if (left_strip_boundary_angle > right_boundary_roi_base_peak ||
            right_strip_boundary_angle < left_boundary_roi_base_peak) {
            continue; // skip strip if not in range
        }

        // maybe calculate once for all modules - or at least save - easier for
        // visualization, debugging but more loops
        double flatfield_normalized_photon_counts =
            (frame.photon_counts(strip_index) + 1) *
            flat_field->get_inverse_normalized_flatfield()(
                strip_index); // we add one to the photon counts to skew
                              // the distribution //TODO obviously
                              // calculate inverse flaftield beforehand

        double some_flatfield_error = 1.0; // TODO: some dummy value - implement

        // I guess it measures the
        // expcected noise - where is the formula
        double photon_counts_variance =
            1. / (std::pow(flatfield_normalized_photon_counts, 2) *
                  (1. / (frame.photon_counts(strip_index) + 1.) +
                   std::pow(some_flatfield_error *
                                flat_field->get_inverse_normalized_flatfield()(
                                    strip_index),
                            2)));

        double strip_width_angle = angular_strip_width_from_DG_parameters(
            DGparameters.centers(module_index),
            DGparameters.conversions(module_index),
            DGparameters.offsets(module_index), local_strip_index);

        double bin_coverage_of_strip = histogram_bin_width / strip_width_angle;

        size_t left_bin_index_covered_by_strip =
            std::max(base_peak_as_bin_index - 50,
                     static_cast<size_t>(
                         (left_strip_boundary_angle + frame.detector_angle) /
                         histogram_bin_width));

        size_t right_bin_index_covered_by_strip =
            std::min(base_peak_as_bin_index + 50,
                     static_cast<size_t>(
                         (right_strip_boundary_angle + frame.detector_angle) /
                         histogram_bin_width));

        // TODO: do strips overlap? - if not second loop is not needed
        for (size_t bin_index = left_bin_index_covered_by_strip;
             bin_index <= right_bin_index_covered_by_strip; ++bin_index) {

            // TODO: ANtonio recalculates bin width and checks if its in
            // boundaries - i think its redundant - check

            // well still weird doesnt change for bins - maybe bin_coverage is
            // wrong
            if (bin_coverage_of_strip >= 0.0001) {

                new_fixed_angle_width_bins_photon_variance(bin_index) +=
                    photon_counts_variance /
                    (bin_coverage_of_strip * bin_coverage_of_strip);

                new_fixed_angle_width_bins_photon_counts(bin_index) +=
                    flatfield_normalized_photon_counts / bin_coverage_of_strip;

                S0(bin_index) +=
                    new_fixed_angle_width_bins_photon_variance(bin_index);
                S1(bin_index) +=
                    new_fixed_angle_width_bins_photon_counts(bin_index) *
                    new_fixed_angle_width_bins_photon_counts(bin_index);
                S2(bin_index) +=
                    new_fixed_angle_width_bins_photon_counts(bin_index) *
                    new_fixed_angle_width_bins_photon_counts(bin_index) *
                    new_fixed_angle_width_bins_photon_counts(bin_index);
            }
        }
    }

    return std::make_tuple(new_fixed_angle_width_bins_photon_counts,
                           new_fixed_angle_width_bins_photon_variance);

    // mmh doesnt it mix up goodness-of-fit criterion
    /*
    for (size_t bin_index = 0; bin_index < num_bins; ++bin_index) {

        new_fixed_angle_width_bin_histogram(bin_index) =
            (new_statistical_weights(bin_index) <=
             std::numeric_limits<double>::epsilon())
                ? 0.
                : new_photon_counts(bin_index) /
                      new_statistical_weights(bin_index);

        // TODO: still need gbinne - maybe its not used
    }
    */
}

// might be deprecated
void AngleCalibration::calculate_fixed_bin_angle_width_histogram(
    const std::vector<std::string> &file_list) {

    new_photon_counts = NDArray<double, 1>(std::array<ssize_t, 1>{num_bins});

    new_photon_count_errors =
        NDArray<double, 1>(std::array<ssize_t, 1>{num_bins});

    // TODO: maybe group these into a 2d array - better cache usage
    NDArray<double, 1> bin_counts(std::array<ssize_t, 1>{num_bins}, 0.0);
    NDArray<double, 1> new_statistical_weights(std::array<ssize_t, 1>{num_bins},
                                               0.0);
    NDArray<double, 1> new_errors(std::array<ssize_t, 1>{num_bins}, 0.0);

    NDView<double, 1> inverse_normalized_flatfield =
        flat_field->get_inverse_normalized_flatfield();

    for (const auto &file : file_list) {
        MythenFrame frame = mythen_file_reader->read_frame(file);
        redistribute_photon_counts_to_fixed_angle_bins(
            frame, bin_counts.view(), new_statistical_weights.view(),
            new_errors.view(), inverse_normalized_flatfield);
    }

    for (ssize_t i = 0; i < new_photon_counts.size(); ++i) {
        new_photon_counts[i] = (new_statistical_weights[i] <=
                                std::numeric_limits<double>::epsilon())
                                   ? 0.
                                   : bin_counts[i] / new_statistical_weights[i];
        new_photon_count_errors[i] =
            (bin_counts[i] <= std::numeric_limits<double>::epsilon())
                ? 0.
                : 1.0 / std::sqrt(bin_counts[i]);
    }
}

NDArray<double, 1> AngleCalibration::calculate_fixed_bin_angle_width_histogram(
    const std::string &file_name) {

    // TODO: maybe group these into a 2d array - better cache usage
    NDArray<double, 1> bin_counts(std::array<ssize_t, 1>{num_bins}, 0.0);
    NDArray<double, 1> new_statistical_weights(std::array<ssize_t, 1>{num_bins},
                                               0.0);
    NDArray<double, 1> new_errors(std::array<ssize_t, 1>{num_bins}, 0.0);

    NDView<double, 1> inverse_normalized_flatfield =
        flat_field->get_inverse_normalized_flatfield();

    MythenFrame frame = mythen_file_reader->read_frame(file_name);

    redistribute_photon_counts_to_fixed_angle_bins(
        frame, bin_counts.view(), new_statistical_weights.view(),
        new_errors.view(), inverse_normalized_flatfield);

    for (ssize_t i = 0; i < bin_counts.size(); ++i) {
        bin_counts[i] = (new_statistical_weights[i] <=
                         std::numeric_limits<double>::epsilon())
                            ? 0.
                            : bin_counts[i] / new_statistical_weights[i];
    }

    return bin_counts;
}

void AngleCalibration::redistribute_photon_counts_to_fixed_angle_bins(
    const MythenFrame &frame, NDView<double, 1> bin_counts,
    NDView<double, 1> new_statistical_weights, NDView<double, 1> new_errors,
    NDView<double, 1> inverse_normalized_flatfield) const {

    ssize_t channel = 0; // TODO handle mask - FlatField still 1d

    if (frame.photon_counts.shape()[0] != mythen_detector->num_strips()) {
        throw std::runtime_error("wrong number of strips read");
    }

    ssize_t num_bins1 = mythen_detector->min_angle() / histogram_bin_width;
    ssize_t num_bins2 = mythen_detector->max_angle() / histogram_bin_width;

    LOG(TLogLevel::logINFO)
        << "position: " << frame.detector_angle << std::endl;

    double exposure_rate = 1. / mythen_detector->exposure_time();

    for (ssize_t strip_index = 0; strip_index < mythen_detector->num_strips();
         ++strip_index) {

        size_t module_index =
            strip_index / MythenDetectorSpecifications::strips_per_module();

        if (mythen_detector->get_bad_channels()[strip_index]) {
            continue;
        }

        double poisson_error = std::sqrt(frame.photon_counts(strip_index)) *
                               inverse_normalized_flatfield(strip_index) *
                               exposure_rate;

        double corrected_photon_count =
            frame.photon_counts(strip_index) *
            inverse_normalized_flatfield(strip_index) * exposure_rate;

        size_t local_strip_index = global_to_local_strip_index_conversion(
            strip_index); // strip_index relative to module

        double diffraction_angle = diffraction_angle_from_DG_parameters(
            DGparameters.centers(module_index),
            DGparameters.conversions(module_index),
            DGparameters.offsets(module_index), local_strip_index);

        diffraction_angle += (frame.detector_angle + mythen_detector->dtt0() +
                              mythen_detector->bloffset());

        if (diffraction_angle < mythen_detector->min_angle() ||
            diffraction_angle > mythen_detector->max_angle())
            continue;

        double angle_covered_by_strip = angular_strip_width_from_DG_parameters(
            DGparameters.centers(module_index),
            DGparameters.conversions(module_index),
            DGparameters.offsets(module_index), local_strip_index);

        double photon_count_per_bin = histogram_bin_width *
                                      corrected_photon_count /
                                      angle_covered_by_strip;
        double error_photon_count_per_bin =
            histogram_bin_width * poisson_error / angle_covered_by_strip;

        double statistical_weights =
            1.0 / std::pow(error_photon_count_per_bin, 2); // 1./sigma²

        double strip_boundary_left =
            diffraction_angle - 0.5 * angle_covered_by_strip;
        double strip_boundary_right =
            diffraction_angle + 0.5 * angle_covered_by_strip;

        ssize_t left_bin_index = std::max(
            num_bins1,
            static_cast<ssize_t>(
                std::floor(strip_boundary_left / histogram_bin_width) - 1));
        ssize_t right_bin_index = std::min(
            num_bins2,
            static_cast<ssize_t>(
                std::ceil(strip_boundary_right / histogram_bin_width) + 1));

        // TODO should it be < or <=
        for (ssize_t bin = left_bin_index; bin <= right_bin_index; ++bin) {
            double bin_coverage = std::min(strip_boundary_right,
                                           (bin + 0.5) * histogram_bin_width) -
                                  std::max(strip_boundary_left,
                                           (bin - 0.5) * histogram_bin_width);

            double bin_coverage_factor = bin_coverage / histogram_bin_width;

            ssize_t bin_index = bin - num_bins1;
            // TODO: maybe have this threshold configurable - or should it
            // be std::numeric_limits
            if (bin_coverage >= 0.0001) {
                new_statistical_weights(bin_index) +=
                    statistical_weights * bin_coverage_factor;
                bin_counts(bin_index) += statistical_weights *
                                         bin_coverage_factor *
                                         photon_count_per_bin;
                new_errors(bin_index) += statistical_weights *
                                         bin_coverage_factor *
                                         std::pow(photon_count_per_bin, 2);
            }
        }
    }
}

// actually used to optimize BC parameters, first parameters used to optimze Lm
// second parameter used to optimize phi
void AngleCalibration::optimization_algorithm(const size_t module_index,
                                              const double shift_parameter1,
                                              const double shift_parameter2) {

    constexpr double tolerance1 =
        0.001; // dont know if this should be configurable

    std::vector<std::pair<double, double>> shift_parameters;
    shift_parameters.reserve(9);

    // TODO dont know if order matters for faster convergence - kept it like in
    // Antonios code
    shift_parameters.emplace_back(-shift_parameter1, -shift_parameter2);
    shift_parameters.emplace_back(-shift_parameter1, shift_parameter2);
    shift_parameters.emplace_back(shift_parameter1, -shift_parameter2);
    shift_parameters.emplace_back(shift_parameter1, shift_parameter2);
    shift_parameters.emplace_back(-shift_parameter1, 0.0);
    shift_parameters.emplace_back(shift_parameter1, 0.0);
    shift_parameters.emplace_back(0.0, -shift_parameter2);
    shift_parameters.emplace_back(0.0, shift_parameter2);
    shift_parameters.emplace_back(0.0, 0.0);

    double previous_similarity_of_peaks =
        calculate_similarity_of_peaks(module_index);

    double next_similarity_of_peaks{};

    std::vector<double> sp(
        shift_parameters
            .size()); // store values of similarity of peak for different shift
                      // parameters //0. sp_x-1,y-1 / 1. sp_x-1,y+1 / 2.
                      // sp_x+1,y-1 / 3. sp_x+1,y+1 / 4. sp_x-1,y / 5. sp_x+1,y
                      // / 6. sp_x,y-1 / 7. sp_x,y+1 / 8 sp_x,y

    bool convergence_criterion = false;
    while (!convergence_criterion) {

        // TODO can i get rid of the break or does it need to be in that order
        // due to convergence
        for (auto shift_parameter : shift_parameters) {
            DGparameters.centers(module_index) += shift_parameter.first;
            DGparameters.conversions(module_index) += shift_parameter.second;
            // TODO: pass parameters directly
            next_similarity_of_peaks =
                calculate_similarity_of_peaks(module_index);
            if (next_similarity_of_peaks < previous_similarity_of_peaks) {
                // update centers for real
                previous_similarity_of_peaks = next_similarity_of_peaks;
                break;
            } else {
                DGparameters.centers(module_index) -= shift_parameter.first;
                DGparameters.conversions(module_index) -=
                    shift_parameter.second;
            }
        }

        // deviates from Antonios code
        // calculate Hessian matrix
        double Dy =
            (sp[7] - sp[6] + sp[1] - sp[0] + sp[3] - sp[2]) /
            (6 * shift_parameter2); //(sp_x,y+1 - sp_x,y-1)/delta_y +
                                    //(sp_x-1,y+1 - sp_x-1,y-1)/delta_y +
                                    //(sp_x+1,y+1 - sp_x+1,y-1)/delta_y
                                    ////averaged central differences
        double Dx =
            (sp[5] - sp[4] + sp[2] - sp[0] + sp[3] - sp[1]) /
            (6 * shift_parameter1); //(sp_x+1,0 - sp_x-1,0)/delta_x +
                                    //(sp_x+1,y-1 - sp_x-1,y-1)/delta_x +
                                    //(sp_x+1,y+1 - sp_x-1, y+1)/2delta_x
        double Dxx =
            ((sp[3] - 2 * sp[7] + sp[1]) + (sp[5] - 2 * sp[8] + sp[4]) +
             (sp[2] - 2 * sp[6] + sp[0])) /
            (shift_parameter1 * shift_parameter2 *
             3); // (sp_x+1,y+1 - 2sp_x,y+1 + sp_x-1,y-1)/delta_x*delta_x etc.
        double Dyy =
            ((sp[3] - 2 * sp[5] + sp[2]) + (sp[7] - 2 * sp[8] + sp[6]) +
             (sp[1] - 2 * sp[4] + sp[0])) /
            (shift_parameter1 * shift_parameter2 * 3);

        double Dyx =
            ((sp[3] - sp[1]) - sp[2] + sp[0]) /
            (4 * shift_parameter1 *
             shift_parameter2); // ((sp_x+1,y+1 - sp_x-1,y+1)/(2*delta_x) -
                                // (sp_x+1,y-1 -
                                // sp_x-1,y-1)/(2*delta_x))/2delta_y

        // calculate eigenvalues of Hessian for regularized inversed Hessian
        double eigenvalue1 =
            0.5 * (Dxx + Dyy) - std::sqrt(0.25 * (Dxx - Dyy) * (Dxx - Dyy) +
                                          Dyx * Dyx); // TODO might be negative?
        double eigenvalue2 =
            0.5 * (Dxx + Dyy) +
            std::sqrt(0.25 * (Dxx - Dyy) * (Dxx - Dyy) + Dyx * Dyx);

        double regularization_term = std::max(
            0.0,
            tolerance1 -
                std::min(eigenvalue1,
                         eigenvalue2)); // Why is this done? what if eigenvalue
                                        // is bigger than tolerance and thus
                                        // value becomes negative?

        double inverse_determinant =
            1. / (Dxx * Dyy - Dyx * Dyx - regularization_term * (Dxx + Dyy) +
                  regularization_term * regularization_term);

        // steepest descent in direction of Hessian^(-1)*gradient Hessian^(-1) =
        // 1/determinant [[Dyy + regularization, -Dyx], [-Dyx, Dxx +
        // regularization]]
        std::pair<double, double> steepest_descent(
            ((Dyy + regularization_term) * Dx - Dyx * Dy) * inverse_determinant,
            ((Dxx + regularization_term) * Dy - Dyx * Dx) *
                inverse_determinant);

        double scale_factor = std::min(
            1.0,
            1. / std::max(std::abs(steepest_descent.first /
                                   DGparameters.centers(module_index)),
                          std::abs(steepest_descent.second /
                                   DGparameters.conversions(module_index))));

        steepest_descent.first *= scale_factor;
        steepest_descent.second *= scale_factor;

        DGparameters.centers(module_index) += steepest_descent.first;
        DGparameters.conversions(module_index) += steepest_descent.second;

        bool found_best_parameters = false;

        while (!found_best_parameters) {

            next_similarity_of_peaks =
                calculate_similarity_of_peaks(module_index);

            found_best_parameters =
                next_similarity_of_peaks <= previous_similarity_of_peaks;
            steepest_descent.first *= 0.5;
            steepest_descent.second *= 0.5;
            DGparameters.centers(module_index) += steepest_descent.first;
            DGparameters.conversions(module_index) += steepest_descent.second;
        }

        double relative_change =
            (previous_similarity_of_peaks - next_similarity_of_peaks) /
            previous_similarity_of_peaks;

        double some_other_criterion = std::sqrt(
            std::pow(Dx * steepest_descent.first + Dy * steepest_descent.second,
                     2) /
            std::pow(steepest_descent.first + steepest_descent.second, 2));

        convergence_criterion =
            (relative_change < 0.001) &
            (some_other_criterion <
             0.001); // TODO: should these tolerances also be configurable -
                     // maybe pass as function argument
    }
}

void AngleCalibration::write_to_file(
    const std::string &filename, const bool store_nonzero_bins,
    const std::filesystem::path &filepath) const {

    std::ofstream output_file(filepath / filename);

    if (!output_file) {
        LOG(TLogLevel::logERROR) << "Error opening file!" << std::endl;
    }

    output_file << std::fixed << std::setprecision(15);

    for (ssize_t i = 0; i < num_bins; ++i) {
        if (new_photon_counts[i] <= std::numeric_limits<double>::epsilon() &&
            store_nonzero_bins) {
            continue;
        }

        output_file << std::floor(mythen_detector->min_angle() /
                                  histogram_bin_width) *
                               histogram_bin_width +
                           i * histogram_bin_width
                    << " " << new_photon_counts[i] << " "
                    << new_photon_count_errors[i] << std::endl;
    }
    output_file.close();
}

} // namespace angcal
