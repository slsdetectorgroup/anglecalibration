/**
 *@brief abstract base classes for simple file interfaces
 * should simplify reading different file formats for angular calibration
 */
#pragma once
#include "aare/Dtype.hpp"
#include "aare/Frame.hpp"
#include "aare/logger.hpp"

/**
 * @brief very simple FileInterface for files that store photon counts of
 * acquisitions
 * @note all functions are pure virtual and must be implemented by the derived
 * classes
 */
class DetectorFileInterface {
  public:
    /**
     * @brief  one frame from the file at the current position
     * @return Frame
     */
    virtual aare::Frame read_frame() = 0;

    /**
     * @brief read one frame from the file at the current position and store it
     * in the provided buffer
     * @param image_buf buffer to store the frame
     * @return void
     */
    virtual void read_into(std::byte *image_buf) = 0;

    virtual void open(const std::string &filename) {
        m_file.close();
        m_filename = filename;

        if (m_mode == "r") {
            try {
                m_file.open(m_filename, std::ios_base::in);
                if (!m_file.good()) {
                    throw std::logic_error("file does not exist");
                }
            } catch (std::exception &e) {

                std::cerr << "Error: " << e.what()
                          << std::endl; // TODO replace with log
            }
        } else {
            throw std::runtime_error(
                LOCATION + "Unsupported mode. Can only CustomMythenFile.");
        }
    }

    /**
     * @brief get the size of one frame in bytes
     * @return size of one frame
     */
    virtual size_t bytes_per_frame() = 0;

    /**
     * @brief get the number of pixels in one frame
     * @return number of pixels in one frame
     */
    virtual size_t pixels_per_frame() = 0;

    /**
     * @brief get the number of rows in the file
     * @return number of rows in the file
     */
    virtual size_t rows() const { return m_rows; }
    /**
     * @brief get the number of columns in the file
     * @return number of columns in the file
     */
    virtual size_t cols() const { return m_cols; };

    virtual ~DetectorFileInterface() { m_file.close(); }

  protected:
    std::string m_mode{};
    ssize_t m_rows{};
    ssize_t m_cols{};
    std::string m_filename{};
    std::ifstream m_file{};
};

/**
 * @brief very simple FileInterface to read from simple ascii files into an
 * array
 * @note all functions are pure virtual and must be implemented by the derived
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
            LOG(aare::TLogLevel::logERROR) << "file does not exist\n";
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