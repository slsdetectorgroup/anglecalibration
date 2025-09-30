#include "AngleCalibration.hpp"

#include "logger.hpp"

#include <csignal>
#include <cstdlib>
#include <filesystem>

#ifdef ANGCAL_PLOT
#include "PlotHelpers.hpp"
#include "plot_histogram.hpp"
#endif

namespace angcal {

inline void handle_sigint(int) {
    LOG(TLogLevel::logINFO) << "\nSIGINT caught — shutting down...\n";
    std::exit(EXIT_SUCCESS);
}
struct SignalHandler {
    SignalHandler() { std::signal(SIGINT, handle_sigint); }
};

AngleCalibration::AngleCalibration(
    std::shared_ptr<MythenDetectorSpecifications> mythen_detector_,
    std::shared_ptr<FlatField> flat_field_,
    std::optional<std::shared_ptr<MythenFileReader>> mythen_file_reader_,
    std::optional<std::shared_ptr<SimpleFileInterface>> custom_file_ptr_)
    : mythen_detector(mythen_detector_), flat_field(flat_field_) {

    if (mythen_file_reader_.has_value()) {
        mythen_file_reader = mythen_file_reader_.value();
    } else {
        mythen_file_reader = std::make_shared<MythenFileReader>();
    }

    if (custom_file_ptr_.has_value()) {
        custom_file_ptr = custom_file_ptr_.value();
    }

    DGparameters = DGParameters(mythen_detector->max_modules());
    BCparameters = BCParameters(mythen_detector->max_modules());
}

void AngleCalibration::set_histogram_bin_width(double bin_width) {
    histogram_bin_width = bin_width;
}

double AngleCalibration::get_histogram_bin_width() const {
    return histogram_bin_width;
}

ssize_t AngleCalibration::get_base_peak_ROI_num_bins() const {
    return 2 * static_cast<ssize_t>(base_peak_roi / histogram_bin_width) + 1;
}

void AngleCalibration::set_base_peak_ROI(const double base_peak_roi_) {
    base_peak_roi = base_peak_roi_;
}

double AngleCalibration::get_base_peak_ROI() const { return base_peak_roi; }

std::shared_ptr<MythenDetectorSpecifications>
AngleCalibration::get_detector_specifications() const {
    return mythen_detector;
}

ssize_t AngleCalibration::new_number_of_bins() const {
    ssize_t new_num_bins =
        std::floor(mythen_detector->max_angle() / histogram_bin_width) -
        std::floor(mythen_detector->min_angle() / histogram_bin_width) + 1;
    return new_num_bins;
}

const DGParameters &AngleCalibration::get_DGparameters() const {
    return DGparameters;
}

const BCParameters &AngleCalibration::get_BCparameters() const {
    return BCparameters;
}

// TODO: command is connected - maybe use that
bool AngleCalibration::module_is_disconnected(const size_t module_index) const {

    // all channels in the module are bad
    return std::all_of(mythen_detector->get_bad_channels().begin() +
                           module_index * mythen_detector->strips_per_module(),
                       mythen_detector->get_bad_channels().begin() +
                           (module_index + 1) *
                               mythen_detector->strips_per_module(),
                       [](const auto &elem) { return elem; });
}

void AngleCalibration::read_initial_calibration_from_file(
    const std::string &filename) {

    custom_file_ptr->open(filename);
    custom_file_ptr->read_into(DGparameters.parameters.buffer(), 8);

    DGparameters.convert_to_BCParameters(
        BCparameters); // initialize BC parameters
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

// TODO: mmh maybe template these on parameter type - use case distinction
double AngleCalibration::diffraction_angle_from_DG_parameters(
    const size_t module_index, const double detector_angle, size_t strip_index,
    const double distance_to_strip) const {

    double offset = DGparameters.offsets(module_index);
    double center = DGparameters.centers(module_index);
    double conversion = DGparameters.conversions(module_index);

    strip_index = std::signbit(conversion)
                      ? MythenDetectorSpecifications::strips_per_module() -
                            strip_index - 1
                      : strip_index; // TODO: are the values sored in reserve?

    return offset +
           180.0 / M_PI *
               (center * std::abs(conversion) -
                std::atan((center - (strip_index + distance_to_strip)) *
                          std::abs(conversion))) +
           detector_angle; // + + mythen_detector->dtt0() +
                           // mythen_detector->bloffset();;
}

double AngleCalibration::diffraction_angle_from_BC_parameters(
    const size_t module_index, const double detector_angle, size_t strip_index,
    const double distance_to_strip) const {

    const double angle_module_center_normal =
        BCparameters.angle_center_module_normal(module_index);
    const double distance_center_sample =
        BCparameters.module_center_sample_distances(module_index);
    const double angle_module_center_beam =
        BCparameters.angle_center_beam(module_index);

    strip_index = std::signbit(distance_center_sample)
                      ? MythenDetectorSpecifications::strips_per_module() -
                            strip_index - 1
                      : strip_index; // TODO: are the values sored in reserve?

    return angle_module_center_beam + angle_module_center_normal -
           180.0 / M_PI *
               std::atan(
                   (std::abs(distance_center_sample) *
                        std::sin(M_PI / 180.0 * angle_module_center_normal) +
                    (MythenDetectorSpecifications::strips_per_module() * 0.5 -
                     strip_index - distance_to_strip) *
                        MythenDetectorSpecifications::pitch()) /
                   (std::abs(distance_center_sample) *
                    std::cos(M_PI / 180.0 * angle_module_center_normal))) +
           detector_angle; // + + mythen_detector->dtt0() +
                           // mythen_detector->bloffset();
}

double AngleCalibration::diffraction_angle_from_EE_parameters(
    const double module_center_distance, const double normal_distance,
    const double angle, const double detector_angle, size_t strip_index,
    const double distance_to_strip) const {

    strip_index = std::signbit(normal_distance)
                      ? MythenDetectorSpecifications::strips_per_module() -
                            strip_index - 1
                      : strip_index;

    return angle -
           180.0 / M_PI *
               std::atan((module_center_distance -
                          MythenDetectorSpecifications::pitch() *
                              (strip_index + distance_to_strip)) /
                         std::abs(normal_distance)) +
           detector_angle; //+ mythen_detector->dtt0() +
                           // mythen_detector->bloffset();
}

// TODO maybe template these on parameter type
double AngleCalibration::angular_strip_width_from_DG_parameters(
    const size_t module_index, const size_t local_strip_index) const {

    return std::abs(diffraction_angle_from_DG_parameters(
                        module_index, 0.0, local_strip_index, 0.5) -
                    diffraction_angle_from_DG_parameters(
                        module_index, 0.0, local_strip_index, -0.5));
}

double AngleCalibration::angular_strip_width_from_BC_parameters(
    const size_t module_index, const size_t local_strip_index) const {

    return std::abs(diffraction_angle_from_BC_parameters(
                        module_index, 0.0, local_strip_index, 0.5) -
                    diffraction_angle_from_BC_parameters(
                        module_index, 0.0, local_strip_index, -0.5));
}

double AngleCalibration::angular_strip_width_from_EE_parameters(
    const double module_center_distance, const double normal_distance,
    const double angle, const size_t local_strip_index) const {

    return std::abs(diffraction_angle_from_EE_parameters(
                        module_center_distance, normal_distance, angle, 0.0,
                        local_strip_index, -0.5) -
                    diffraction_angle_from_EE_parameters(
                        module_center_distance, normal_distance, angle, 0.0,
                        local_strip_index, 0.5));
}

void AngleCalibration::set_base_peak_angle(const double base_peak_angle_) {
    base_peak_angle = base_peak_angle_;
}

double AngleCalibration::get_base_peak_angle() const { return base_peak_angle; }

/*
void AngleCalibration::calculate_average_photon_counts() {
    NDArray<size_t, 1> num_frames_with_base_peak_in_module(
        std::array<ssize_t, 1>{mythen_detector->max_modules()});
    NDArray<double, 1> average_photon_counts(
        std::array<ssize_t, 1>{mythen_detector->max_modules()});

    for (const auto &file : file_list) {
        MythenFrame frame = mythen_file_reader->read_frame(file);

        for (size_t module_index; module_index < mythen_detector->max_modules();
             ++module_index) {
            if (base_peak_is_in_module(module_index, frame.detector_angle)) {
                num_frames_with_base_peak_in_module(module_index) += 1;
                average_photon_counts(module_index) = std::accumulate(
                    frame.photon_counts.begin() +
                        module_index * mythen_detector->strips_per_module(),
                    frame.photon_counts.end() +
                        (module_index + 1) *
                            mythen_detector->strips_per_module(),
                    1.0); // mmh why 1.0 though?
            }
        }
    }

    // scale with average over all runs
}
*/

double AngleCalibration::similarity_criterion(const NDView<double, 1> S0,
                                              const NDView<double, 1> S1,
                                              const NDView<double, 1> S2,
                                              const size_t num_runs) const {

    double similarity_criterion = 0;
    for (ssize_t bin_index = 0; bin_index < S0.size(); ++bin_index) {

        double weighted_average =
            S0(bin_index) < std::numeric_limits<double>::epsilon()
                ? 0.0
                : 1. / S0(bin_index); // photon variance over each run //TODO
                                      // check if this is correct with Antonios
                                      // code!! - how to caluclate varaince?
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
           std::max(num_runs,
                    static_cast<size_t>(1)); // should actually only divide by
                                             // runs or also by num bins?
}

// TODO maybe inline this
bool AngleCalibration::base_peak_is_in_module(
    const size_t module_index, const double detector_angle,
    std::optional<double> bounds_in_angles) const {

    if (!bounds_in_angles.has_value()) {
        bounds_in_angles = base_peak_roi; // take ROI of base peak in angles
    }

    double left_module_boundary_angle =
        diffraction_angle_from_DG_parameters(module_index, detector_angle, 0);
    double right_module_boundary_angle = diffraction_angle_from_DG_parameters(
        module_index, detector_angle, mythen_detector->strips_per_module() - 1);

    LOG(TLogLevel::logDEBUG1) << fmt::format(
        "module_boundaries_in_angle for module {} [{}, {}]\n", module_index,
        left_module_boundary_angle, right_module_boundary_angle);

    // TODO check bounds with base_peak_hwid - TODO: I changed it such that it
    // needs to be fully in module - otherwise similarity criterion gets
    // distorted
    return (base_peak_angle + *bounds_in_angles < right_module_boundary_angle &&
            left_module_boundary_angle < base_peak_angle - *bounds_in_angles);
}

double
AngleCalibration::rate_correction_factor(const double photon_count,
                                         const double exposure_time) const {
    const double dead_time = 2.915829802160547e-7; // measured dead-time

    const double maximum_count_rate =
        std::exp(-1); // theoretical maximum count rate for
                      // dead_time*measured_photon_counts_per_second

    double photon_counts_per_second = photon_count / exposure_time;

    photon_counts_per_second =
        std::min(maximum_count_rate, photon_counts_per_second *
                                         dead_time); // multiply width dead time

    // -actual_count_rate*dead_time = W -> -dead_time*measured_count_rate =
    // We^{W} -> W: Lambert W function

    // numerical calculation of inverse of Lambart W - function
    double W_prev_iter = -photon_counts_per_second; // initial guess
    double W_next_iter{};
    bool method_converged = false;
    while (!method_converged) {
        W_next_iter = photon_counts_per_second * std::exp(W_prev_iter);
        method_converged = std::abs(W_next_iter - W_prev_iter) <
                           10 * std::numeric_limits<double>::epsilon();
        W_prev_iter = W_next_iter;
    }

    double actual_count_rate = -W_prev_iter;

    return actual_count_rate /
           photon_counts_per_second; // TODO: its a statistic variable not a
                                     // constant do I take it into account?
}

std::pair<double, double> AngleCalibration::calculate_corrected_photon_counts(
    const double photon_counts, const size_t global_strip_index) const {

    double flatfield_normalized_photon_counts =
        (photon_counts + 1) *
        flat_field->get_inverse_normalized_flatfield()(global_strip_index,0);

    double some_flatfield_error =
        1.0; // TODO: some dummy value - implement read from file

    // I guess it measures the
    // expcected noise - where is the formula - used as the variance
    double inverse_photon_counts_variance =
        1. / (std::pow(flatfield_normalized_photon_counts, 2) *
              (1. / (photon_counts + 1) +
               std::pow(some_flatfield_error *
                            flat_field->get_inverse_normalized_flatfield()(
                                global_strip_index,0),
                        2))); // this is probably the inverse - for easier
                              // multiplication - division by variance

    return std::pair(flatfield_normalized_photon_counts,
                     inverse_photon_counts_variance);
}

double
AngleCalibration::calculate_similarity_of_peaks(const size_t module_index,
                                                PlotHandle plot) const {

    ssize_t num_bins_in_ROI =
        2 * static_cast<ssize_t>(base_peak_roi / histogram_bin_width) + 1;
    // used to calculate similarity criterion between peaks of different
    // acquisition S_index = sum_i^num_runs
    // photon_count^index*photon_variance
    NDArray<double, 1> S2(std::array<ssize_t, 1>{num_bins_in_ROI}, 0.0);
    NDArray<double, 1> S1(std::array<ssize_t, 1>{num_bins_in_ROI}, 0.0);
    NDArray<double, 1> S0(std::array<ssize_t, 1>{num_bins_in_ROI}, 0.0);

#ifdef ANGCAL_PLOT
    auto data_file_path =
        std::filesystem::current_path().parent_path() / "build" / "data";
    if (!std::filesystem::exists(data_file_path))
        std::filesystem::create_directories(data_file_path);
    plot->clear();
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
#endif

    size_t num_runs = 0;
    for (const auto &file : file_list) {

        LOG(TLogLevel::logDEBUG1) << "file: " << file;

        MythenFrame frame = mythen_file_reader->read_frame(file);

        // base peak angle is in module
        if (base_peak_is_in_module(module_index, frame.detector_angle)) {

            LOG(TLogLevel::logDEBUG1) << "file name: " << file;

            NDArray<double, 1> fixed_angle_width_bins_photon_counts(
                std::array<ssize_t, 1>{num_bins_in_ROI}, 0.0);
            NDArray<double, 1>
                inverse_fixed_angle_width_bins_photon_counts_variance(
                    std::array<ssize_t, 1>{num_bins_in_ROI}, 0.0);
            // calculates flatfield normalized photon counts and
            // photon_count_variance for ROI around base_peak and
            // redistributes to fixed angle width bins

            redistribute_photon_counts_to_fixed_angle_width_bins<true>(
                module_index, frame,
                fixed_angle_width_bins_photon_counts.view(),
                inverse_fixed_angle_width_bins_photon_counts_variance.view(),
                S0.view(), S1.view(), S2.view());

#ifdef ANGCAL_PLOT
            auto bin_to_diffraction_angle_base_peak_ROI_only =
                [this](const size_t bin_index) {
                    return bin_index * this->histogram_bin_width -
                           this->base_peak_roi + this->base_peak_angle;
                };

            std::string filename =
                std::filesystem::path(file).stem().string() + ".dat";
            auto dataset_name = data_file_path / filename;

            plot->append_to_plot(
                fixed_angle_width_bins_photon_counts.view(),
                {0,
                 2 * static_cast<ssize_t>(base_peak_roi / histogram_bin_width) +
                     1},
                bin_to_diffraction_angle_base_peak_ROI_only, dataset_name);
            // plot->pause();
#endif
            ++num_runs;
        }
    }

    // handle those cases
    if (num_runs == 0) {
        LOG(TLogLevel::logDEBUG) << "there was no frame where base peak ROI "
                                    "overlapped with module region";
    }

#ifdef ANGCAL_PLOT
    plot->flush(); // make plot appear
    std::this_thread::sleep_for(
        std::chrono::milliseconds(1000)); // let gnuplot update
#endif

    // TODO: give a warning if there was no acquisition for a module

    double similarity_of_peaks =
        num_runs == 0
            ? std::numeric_limits<double>::infinity()
            : similarity_criterion(S0.view(), S1.view(), S2.view(), num_runs);
    LOG(TLogLevel::logDEBUG1) << "similarity_of_peaks: " << similarity_of_peaks;
    return similarity_of_peaks;
}

// TODO: maybe have a function where we calculate the average photon counts
// - multiple loops - or directly store normalized data somewhere instead of
// computing on the fly

void AngleCalibration::calibrate(
    const std::vector<std::string> &file_list_, const double base_peak_angle_,
    std::optional<std::filesystem::path> output_file) {

    file_list = file_list_;

    base_peak_angle = base_peak_angle_;

    std::ofstream file;
    if (output_file.has_value()) {
        file.open(output_file.value(), std::ios::out | std::ios::app);
        if (file.tellp() == 0) {
            file << "module,center,conversion,offset\n"; // only write header if
                                                         // file is empty
        }
    }

    for (size_t module_index = 0; module_index < mythen_detector->max_modules();
         ++module_index) {

        // skip if module is not connected
        if (!module_is_disconnected(module_index)) {

            LOG(angcal::TLogLevel::logINFO)
                << "starting calibration for module " << module_index;

#ifdef ANGCAL_PLOT
            std::string plot_title =
                fmt::format("Base Peaks for module {} ", module_index);
            optimization_algorithm(
                module_index,
                std::make_shared<PlotCalibrationProcess>(plot_title));
#else
            optimization_algorithm(module_index);
#endif
        } else {
            LOG(TLogLevel::logINFO)
                << fmt::format("module {} is disconnected", module_index);
        }

        // write to file
        if (output_file.has_value()) {
            auto [center, conversion, offset] =
                BCparameters.convert_to_DGParameters(module_index);
            append_to_file(file, module_index, center, conversion, offset);
        }
    }

    if (output_file.has_value()) {
        file.close();
    }
}

void AngleCalibration::calibrate(const std::vector<std::string> &file_list_,
                                 const double base_peak_angle_,
                                 const size_t module_index) {

    file_list = file_list_;
    base_peak_angle = base_peak_angle_;

    if (!module_is_disconnected(module_index)) {
        LOG(angcal::TLogLevel::logINFO)
            << "starting calibration for module " << module_index;

#ifdef ANGCAL_PLOT
        std::string plot_title =
            fmt::format("Base Peaks for module {} ", module_index);
        optimization_algorithm(
            module_index, std::make_shared<PlotCalibrationProcess>(plot_title));
#else
        optimization_algorithm(module_index);
#endif
    }
}

NDArray<double, 1>
AngleCalibration::redistribute_photon_counts_to_fixed_angle_width_bins(
    const MythenFrame &frame, const size_t module_index) const {

    ssize_t new_num_bins = new_number_of_bins();

    NDArray<double, 1> fixed_angle_width_bins_photon_counts =
        NDArray<double, 1>(std::array<ssize_t, 1>{new_num_bins}, 0.0);

    NDArray<double, 1> fixed_angle_width_bins_photon_variance =
        NDArray<double, 1>(std::array<ssize_t, 1>{new_num_bins}, 0.0);

    redistribute_photon_counts_to_fixed_angle_width_bins<false>(
        module_index, frame, fixed_angle_width_bins_photon_counts.view(),
        fixed_angle_width_bins_photon_variance.view());

    return fixed_angle_width_bins_photon_counts;
};

NDArray<double, 1>
AngleCalibration::redistribute_photon_counts_to_fixed_angle_width_bins(
    const MythenFrame &frame) const {
    ssize_t new_num_bins = new_number_of_bins();

    NDArray<double, 1> fixed_angle_width_bins_photon_counts =
        NDArray<double, 1>(std::array<ssize_t, 1>{new_num_bins}, 0.0);

    NDArray<double, 1> fixed_angle_width_bins_photon_variance =
        NDArray<double, 1>(std::array<ssize_t, 1>{new_num_bins}, 0.0);

    // TODO: actually they should not be added up each set of modules is
    // independant - at beamline the module positions overlap (e.g.
    // counterclockwise)
    //  - how to handle in a generic way? - depends on detector arrangement
    for (size_t module_index = 0; module_index < mythen_detector->max_modules();
         ++module_index) {

        if (module_is_disconnected(module_index)) {
            continue;
        }

        LOG(TLogLevel::logDEBUG1)
            << fmt::format("module_index {} contibutes", module_index);

        redistribute_photon_counts_to_fixed_angle_width_bins<false>(
            module_index, frame, fixed_angle_width_bins_photon_counts.view(),
            fixed_angle_width_bins_photon_variance.view());
    }

    return fixed_angle_width_bins_photon_counts;
}

NDArray<double, 1>
AngleCalibration::redistributed_photon_counts_in_base_peak_ROI(
    const MythenFrame &frame, const size_t module_index) const {

    NDArray<double, 1> fixed_angle_width_bins_photon_counts =
        NDArray<double, 1>(std::array<ssize_t, 1>{get_base_peak_ROI_num_bins()},
                           0.0);

    NDArray<double, 1> fixed_angle_width_bins_photon_variance =
        NDArray<double, 1>(std::array<ssize_t, 1>{get_base_peak_ROI_num_bins()},
                           0.0);

    redistribute_photon_counts_to_fixed_angle_width_bins<true>(
        module_index, frame, fixed_angle_width_bins_photon_counts.view(),
        fixed_angle_width_bins_photon_variance.view());

    return fixed_angle_width_bins_photon_counts;
}

// actually used to optimize BC parameters, first parameters used to optimze
// Lm second parameter used to optimize phi
void AngleCalibration::optimization_algorithm(const size_t module_index,
                                              PlotHandle gp,
                                              const double shift_parameter1,
                                              const double shift_parameter2) {

    constexpr double tolerance1 =
        0.001; // dont know if this should be configurable

    constexpr size_t max_iterations = 200;

    std::vector<std::pair<double, double>> shift_parameters;
    shift_parameters.reserve(9);

    // TODO dont know if order matters for faster convergence - kept it like
    // in Antonios code
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
        calculate_similarity_of_peaks(module_index, gp);

    double next_similarity_of_peaks{};

    std::vector<double> sp(
        shift_parameters
            .size()); // store values of similarity of peak for different
                      // shift parameters //0. sp_x-1,y-1 / 1. sp_x-1,y+1
                      // / 2. sp_x+1,y-1 / 3. sp_x+1,y+1 / 4. sp_x-1,y / 5.
                      // sp_x+1,y / 6. sp_x,y-1 / 7. sp_x,y+1 / 8 sp_x,y

    bool convergence_criterion = false;
    size_t iteration_index = 0;
    while (!convergence_criterion && (iteration_index < max_iterations)) {

        // TODO can i get rid of the break or does it need to be in that
        // order due to convergence

        for (size_t parameter_index = 0;
             parameter_index < shift_parameters.size(); ++parameter_index) {

            LOG(TLogLevel::logINFO)
                << fmt::format("Iteration {} with peak similarity of {}",
                               iteration_index, previous_similarity_of_peaks);

            ++iteration_index;

            BCparameters.angle_center_module_normal(module_index) +=
                shift_parameters[parameter_index].first;
            BCparameters.module_center_sample_distances(module_index) +=
                shift_parameters[parameter_index].second;
            // TODO: pass parameters directly
            next_similarity_of_peaks =
                calculate_similarity_of_peaks(module_index, gp);
            sp[parameter_index] = next_similarity_of_peaks;
            /*
            BCparameters.angle_center_module_normal(module_index) -=
                shift_parameters[parameter_index].first;
            BCparameters.module_center_sample_distances(module_index) -=
                shift_parameters[parameter_index].second;
            */

            if (next_similarity_of_peaks < previous_similarity_of_peaks) {
                // update centers for real
                previous_similarity_of_peaks = next_similarity_of_peaks;
                // parameter_index = 0;
                //  break;
            } else {
                BCparameters.angle_center_module_normal(module_index) -=
                    shift_parameters[parameter_index].first;
                BCparameters.module_center_sample_distances(module_index) -=
                    shift_parameters[parameter_index].second;
            }
        }

        LOG(TLogLevel::logINFO)
            << "peak similarity: " << previous_similarity_of_peaks;

#ifdef ANGCAL_PLOT
        PlotHelper::pause(); // dont know if handled properly more for debugging
#endif
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
             3); // (sp_x+1,y+1 - 2sp_x,y+1 + sp_x-1,y-1)/delta_x*delta_x
                 // etc.
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

        // steepest descent in direction of Hessian^(-1)*gradient
        // Hessian^(-1) = 1/determinant [[Dyy + regularization, -Dyx],
        // [-Dyx, Dxx + regularization]]
        std::pair<double, double> steepest_descent(
            ((Dyy + regularization_term) * Dx - Dyx * Dy) * inverse_determinant,
            ((Dxx + regularization_term) * Dy - Dyx * Dx) *
                inverse_determinant);

        LOG(TLogLevel::logINFO) << fmt::format(
            " fine tune parameters in direction of steepest descent: [{},{}]",
            steepest_descent.first, steepest_descent.second);

        double scale_factor = std::min(
            1.0,
            1. / std::max(std::abs(steepest_descent.first /
                                   BCparameters.angle_center_module_normal(
                                       module_index)),
                          std::abs(steepest_descent.second /
                                   BCparameters.module_center_sample_distances(
                                       module_index))));

        steepest_descent.first *= scale_factor;
        steepest_descent.second *= scale_factor;

        BCparameters.angle_center_module_normal(module_index) +=
            steepest_descent.first;
        BCparameters.module_center_sample_distances(module_index) +=
            steepest_descent.second;

        bool found_best_parameters = false;

        // check termination criterion iteration_index added on purpose
        while (!found_best_parameters && iteration_index < max_iterations) {

            next_similarity_of_peaks =
                calculate_similarity_of_peaks(module_index, gp);

            // TODO what if im stuck in here
            if (next_similarity_of_peaks > previous_similarity_of_peaks) {
                BCparameters.angle_center_module_normal(module_index) -=
                    steepest_descent.first;
                BCparameters.module_center_sample_distances(module_index) -=
                    steepest_descent.second;
            }

            found_best_parameters =
                next_similarity_of_peaks <= previous_similarity_of_peaks;
            // TODO what if im stuck in here
            steepest_descent.first *= 0.5;
            steepest_descent.second *= 0.5;
            BCparameters.angle_center_module_normal(module_index) +=
                steepest_descent.first;
            BCparameters.module_center_sample_distances(module_index) +=
                steepest_descent.second;

            LOG(TLogLevel::logINFO)
                << fmt::format("Iteration {} with peak similarity of {}",
                               iteration_index, previous_similarity_of_peaks);

            ++iteration_index;
        }

        double relative_change =
            (previous_similarity_of_peaks - next_similarity_of_peaks) /
            previous_similarity_of_peaks;

        LOG(TLogLevel::logINFO)
            << "relative change of similarity of peaks: " << relative_change;

        double some_other_criterion = std::sqrt(
            std::pow(Dx * steepest_descent.first + Dy * steepest_descent.second,
                     2) /
            std::pow(steepest_descent.first + steepest_descent.second, 2));

        convergence_criterion = (relative_change < 0.001); // &
        //(some_other_criterion <
        // 0.001); // TODO: should these tolerances also be configurable -
        // maybe pass as function argument

        ++iteration_index;
    }
}

// TODO: it just writes to a csv/txt file maybe consider using same file format
// as in m_custom_file_ptr but then user needs to implement more e.g. custom
// append
void AngleCalibration::write_to_file(
    const std::filesystem::path &filename) const {

    std::ofstream output_file(filename);

    if (!output_file) {
        LOG(angcal::TLogLevel::logERROR) << "Error opening file!" << std::endl;
    }

    output_file.precision(15);

    output_file << "module, center, conversion, offset" << std::endl; // header

    for (ssize_t module_index = 0; module_index < BCparameters.num_modules();
         ++module_index) {
        auto [center, conversion, offset] =
            BCparameters.convert_to_DGParameters(module_index);
        append_to_file(output_file, module_index, center, conversion, offset);
    }

    output_file.close();
}

void AngleCalibration::append_to_file(std::ofstream &os,
                                      const size_t module_index,
                                      const double center,
                                      const double conversion,
                                      const double offset) const {
    os << module_index << "," << center << "," << conversion << "," << offset
       << std::endl;
}

} // namespace angcal
