#ifndef SCARASOLUTION_H
#define SCARASOLUTION_H
//#include "libs/Module.h"
#include "BaseSolution.h"
#include "Gcode.h"

class Config;

class SCARASolution : public BaseSolution {
    public:
        SCARASolution(Config*);
        void cartesian_to_actuator(const float[], ActuatorCoordinates &) const override;
        void actuator_to_cartesian(const ActuatorCoordinates &, float[] ) const override;

        bool set_optional(const arm_options_t& options) override;
        bool get_optional(arm_options_t& options, bool force_all) const override;

    private:
        void init();
        float to_degrees(float radians) const;
        float to_radians(float degrees) const;

        float arm1_length;
        float arm2_length;
        float scara_offset_x;
        float scara_offset_y;
        float scara_scaling_x;
        float scara_scaling_y;
        float slow_rate;
        bool scara_right_hand;
        bool inverse_psi;
};

#endif // SCARASOLUTION_H
