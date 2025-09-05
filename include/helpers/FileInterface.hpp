/**
 *@brief abstract base classes for simple file interfaces
 * should simplify reading different file formats for angular calibration
 */
#pragma once
#include "aare/Dtype.hpp"
#include "aare/Frame.hpp"
#include "logger.hpp"

/**
 * @brief Simple FileInterface to read from simple ascii files or raw files
 * into a bytestream
 * @note some functions are pure virtual and must be implemented by the derived
 * classes
 */
class SimpleFileInterface {
  public:
    /**
     * @brief read file into an array
     * @param array: view of the array
     * @return void
     */
    virtual void open(const std::string &filename) {
        m_filename = filename;
        m_file.close();
        try {
            m_file = std::ifstream(filename, std::ios_base::in);

            if (!m_file.good()) {
                throw std::logic_error("file does not exist");
            }
        } catch (std::exception &e) {
            LOG(angcal::TLogLevel::logERROR) << "file does not exist\n";
        }
    }

    /**
     * @brief read into byte buffer
     * @param data_type_bytes: if buffer from a vector or array provide number
     * of bytes for datatype
     */
    virtual void read_into(std::byte *image_buf,
                           const ssize_t data_type_bytes = 1) = 0;

    virtual ~SimpleFileInterface() { m_file.close(); }

  protected:
    std::string m_filename{};
    std::ifstream m_file{};
};