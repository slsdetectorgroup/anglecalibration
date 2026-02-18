#pragma once

#include <filesystem>
#include <fstream>

class LogFile {

  public:
    LogFile(std::filesystem::path log_file_path)
        : m_log_file_path(log_file_path) {

        open();
    }

    ~LogFile() {
        if (m_log_file.is_open()) {
            m_log_file.close();
        }
    }

    void close() {
        if (m_log_file.is_open()) {
            m_log_file.close();
        }
    }

    void open() {
        m_log_file.open(m_log_file_path, std::ios::out);

        if (!m_log_file.is_open()) {
            throw std::runtime_error("Could not open log file: " +
                                     m_log_file_path.string());
        }
    }

    void append(const std::string &log_entry) {
        if (m_log_file.is_open()) {
            m_log_file << log_entry;
        } else {
            throw std::runtime_error("Log file is not open: " +
                                     m_log_file_path.string());
        }
    }

  private:
    std::ofstream m_log_file;
    std::filesystem::path m_log_file_path;
};

inline LogFile FlatField_correctionfactor("FlatField_correctionfactor.log");
inline LogFile IncidentIntensity_correctionfactor(
    "IncidentIntensity_correctionfactor.log");
inline LogFile Rate_correctionfactor("Rate_correctionfactor.log");
inline LogFile Transversewidth_correction("Transversewidth_correction.log");

inline LogFile FlatFieldErrors("FlatFieldErrors.log");

inline LogFile CorrectedPhotonCountsLogFile("CorrectedPhotonCounts.log");

inline LogFile
    CorrectedPhotonCountsErrorsLogFile("CorrectedPhotonCountsErrors.log");

inline LogFile BInBoundariesLogFile("BinBoundaries.log");

inline LogFile StripAngles("StripAngles.log");

inline LogFile CenterStripAngles("CenterStripAngles.log");

inline LogFile EnocderMs("EncoderMs.log");