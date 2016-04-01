//
// Created by Michael Afanasiev on 2016-03-16.
//

#ifndef SALVUS_FIELDIO_H
#define SALVUS_FIELDIO_H

#include <Eigen/Dense>
#include <H5public.h>
#include <H5Ipublic.h>
#include "../Utilities/Options.h"

/**
 * Base class which implements IO operations for things like strain, velocity, etc.
 */
class FieldIo {

protected:

    int mPolynomialDegree;
    std::string mPhysicsSystem;

public:

    /**
     * Base constructor setting generic options (polynomial degree, physics).
     * @param [in] options User defined options class.
     */
    FieldIo(Options options);

    /**
     * Virtual desctructor. Meant to be overridden by the derived classes.
     */
    virtual ~FieldIo() {};

    /**
     * Factory which returns a different instance of FieldIo, depending on
     * user defined options.
     */
    static FieldIo *factory(Options options);

    /**
     * Saves a field variable at a given time.
     * @param [in] val Eigen::MatrixXd containing the field to save.
     * @param [in] time Simulation time.
     */
    virtual void save(const Eigen::MatrixXd &val, const double &time) = 0;

    /**
     * Reads in a field variable at a given time.
     * @param[in] time Simulation time.
     */
    virtual void read(const double &time) = 0;

    /**
     * Takes any buffered saves and writes them to disk.
     * @param [in] time Simulation time.
     */
    virtual void write(const double &time) = 0;

};

class FieldIoHdf5 : public FieldIo {

private:

    struct snapshot {
        std::vector<float>      time;
        std::vector<char>       nByte;
        std::vector<char>       byte1;
        std::vector<short int>  byte2;
        std::vector<float>      byte4;
        std::vector<double>     byte8;
        std::vector<float>      scale;
        std::vector<float>      offset;
    };
    snapshot mSnapshot;

    hid_t mPlistId, mFileId;
    size_t mTsteps;

    void _write_disk();

public:

    FieldIoHdf5(Options options);
    ~FieldIoHdf5();
    virtual void save(const Eigen::MatrixXd &val, const double &time);
    virtual void read(const double &time);
    virtual void write(const double &time);

};


#endif //SALVUS_FIELDIO_H
