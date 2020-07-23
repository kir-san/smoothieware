#include "SCARASolution.h"
#include "ActuatorCoordinates.h"
#include "ConfigValue.h"
#include "StreamOutputPool.h"
#include "checksumm.h"
#include "libs/Kernel.h"
#include <fastmath.h>

#include "libs/nuts_bolts.h"

#include "libs/Config.h"

#define arm1_length_checksum CHECKSUM("arm1_length")
#define arm2_length_checksum CHECKSUM("arm2_length")
#define scara_offset_x_checksum CHECKSUM("scara_offset_x")
#define scara_offset_y_checksum CHECKSUM("scara_offset_y")
#define scara_scaling_x_checksum CHECKSUM("scara_scaling_x")
#define scara_scaling_y_checksum CHECKSUM("scara_scaling_y")
#define scara_homing_checksum CHECKSUM("scara_homing")
#define scara_right_hand_checksum CHECKSUM("scara_right_hand")
#define inverse_psi_checksum CHECKSUM("inverse_psi")

#define SQ(x) powf(x, 2)
#define ROUND(x, y) (roundf(x * 1e##y) / 1e##y)

SCARASolution::SCARASolution(Config *config)
{
    // arm1_length is the length of the inner main arm from hinge to hinge
    arm1_length = config->value(arm1_length_checksum)->by_default(150.0f)->as_number();
    // arm2_length is the length of the inner main arm from hinge to hinge
    arm2_length = config->value(arm2_length_checksum)->by_default(150.0f)->as_number();
    // scara_offset_x is the x offset of bed zero position towards the SCARA tower center
    scara_offset_x = config->value(scara_offset_x_checksum)->by_default(100.0f)->as_number();
    // scara_offset_y is the y offset of bed zero position towards the SCARA tower center
    scara_offset_y = config->value(scara_offset_y_checksum)->by_default(-60.0f)->as_number();
    // Axis scaling is used in final calibration
    scara_scaling_x = config->value(scara_scaling_x_checksum)->by_default(1.0F)->as_number(); // 1 = 100% : No scaling
    scara_scaling_y = config->value(scara_scaling_y_checksum)->by_default(1.0F)->as_number();
    // scara_undefined is the ratio at which the SCARA position is undefined.
    // required to prevent the arm moving through singularity points
    scara_right_hand = config->value(scara_right_hand_checksum)->by_default(true)->as_bool();

    init();
}

void SCARASolution::init()
{
}

float SCARASolution::to_degrees(float radians) const
{
    return radians * (180.0F / 3.14159265359f);
}

float SCARASolution::to_radians(float radians) const
{
    return radians / (180.0F / 3.14159265359f);
}

void SCARASolution::cartesian_to_actuator(const float cartesian_mm[], ActuatorCoordinates &actuator_mm) const
{
    // THEKERNEL->streams->printf("DEBUG: cartesian_to_actuator (start)\n");

    float
        SCARA_pos[2], // Хранение координат
        SCARA_C2;     // Переменная позволяющая определить способ расчета углов

    // Конвертирование текущих координат в координаты Скары, учитывая смещение башни и множитель размера
    SCARA_pos[X_AXIS] = cartesian_mm[X_AXIS] * this->scara_scaling_x - this->scara_offset_x;
    SCARA_pos[Y_AXIS] = (cartesian_mm[Y_AXIS] * this->scara_scaling_y - this->scara_offset_y);

    // Если плечи одной длины, то используется одна формула расчета, если разной, то другая
    if (this->arm1_length == this->arm2_length)
        SCARA_C2 = (SQ(SCARA_pos[X_AXIS]) + SQ(SCARA_pos[Y_AXIS])) / (2.0f * SQ(this->arm1_length)) - 1;
    else
        SCARA_C2 = (SQ(SCARA_pos[X_AXIS]) + SQ(SCARA_pos[Y_AXIS]) - SQ(this->arm1_length) - SQ(this->arm2_length)) / (2.0f * this->arm1_length * this->arm2_length);

    // При выходе за пределы коэффициентов, возникают ошибки в расчетах и результаты не соответствуют действительности
    // Если значение переменной вышло за пределы коэффициентов, то используется тригонометрия
    if (SCARA_C2 > 0.95f || SCARA_C2 < -0.95f)
    {
        float psi = acos(SCARA_C2);

        actuator_mm[BETA_STEPPER] = to_degrees(psi); // Угол для Y

        float arctan_xy = atan2f(SCARA_pos[Y_AXIS], SCARA_pos[X_AXIS]);
        float arctan = atan2f(this->arm2_length * sin(psi),
                              this->arm1_length + this->arm2_length * cos(psi));

        actuator_mm[ALPHA_STEPPER] = to_degrees(arctan_xy - arctan); // Угол для X
    }
    // Если не вышло, то используется алгебра так как быстрее
    else
    {
        float s2 = sqrtf(1.0f - SQ(SCARA_C2));
        float k1 = this->arm1_length + this->arm2_length * SCARA_C2;
        float k2 = this->arm2_length * s2;
        float theta = (atan2f(SCARA_pos[X_AXIS], SCARA_pos[Y_AXIS]) - atan2f(k1, k2)) * -1.0f;
        actuator_mm[ALPHA_STEPPER] = to_degrees(theta);
        float psi = atan2f(s2, SCARA_C2);
        actuator_mm[BETA_STEPPER] = to_degrees(psi);
    }

    // THEKERNEL->streams->printf("DEBUG: cartesian_mm x=%f, y=%f, z=%f\n", cartesian_mm[X_AXIS], cartesian_mm[Y_AXIS], cartesian_mm[Z_AXIS]);
    // THEKERNEL->streams->printf("DEBUG: SCARA_pos x=%f, y=%f\n", SCARA_pos[X_AXIS], SCARA_pos[Y_AXIS]);

    // Если скара леворукая, то все углы меняют знак
    if (!scara_right_hand)
    {
        actuator_mm[ALPHA_STEPPER] = actuator_mm[ALPHA_STEPPER] * (-1.0f);
        actuator_mm[BETA_STEPPER] = actuator_mm[BETA_STEPPER] * (-1.0f);
    }

    actuator_mm[GAMMA_STEPPER] = cartesian_mm[Z_AXIS]; // No inverse kinematics on Z - Position to add bed offset?

    // THEKERNEL->streams->printf("DEBUG: actuator_mm x=%f, y=%f, z=%f\n", actuator_mm[ALPHA_STEPPER], actuator_mm[BETA_STEPPER], actuator_mm[GAMMA_STEPPER]);
    // THEKERNEL->streams->printf("DEBUG: cartesian_to_actuator (end)\n");
}

void SCARASolution::actuator_to_cartesian(const ActuatorCoordinates &actuator_mm, float cartesian_mm[]) const
{
    // Perform forward kinematics, and place results in cartesian_mm[]
    float x_sin, x_cos, y_sin, y_cos, actuator_rad[2];

    actuator_rad[X_AXIS] = to_radians(actuator_mm[X_AXIS]);
    actuator_rad[Y_AXIS] = to_radians(actuator_mm[Y_AXIS] + actuator_mm[X_AXIS]);

    if (inverse_psi)
        actuator_rad[Y_AXIS] = to_radians(180.0f - actuator_mm[Y_AXIS] + actuator_mm[X_AXIS]);
    ;

    x_sin = sinf(actuator_rad[X_AXIS]) * this->arm1_length;
    x_cos = cosf(actuator_rad[X_AXIS]) * this->arm1_length;

    y_sin = sinf(actuator_rad[Y_AXIS]) * this->arm2_length;
    y_cos = cosf(actuator_rad[Y_AXIS]) * this->arm2_length;

    cartesian_mm[X_AXIS] = (x_cos + y_cos + scara_offset_x) / scara_scaling_x; // Theta
    if (scara_right_hand)
    {
        cartesian_mm[Y_AXIS] = (x_sin + y_sin + scara_offset_y) / scara_scaling_y; //psi
    }
    else
    {
        cartesian_mm[Y_AXIS] = (scara_offset_y - x_sin - y_sin) / scara_scaling_y; //psi
    }

    cartesian_mm[Z_AXIS] = actuator_mm[Z_AXIS];

    // THEKERNEL->streams->printf("DEBUG: cartesian_mm[X_AXIS] = %f\n", cartesian_mm[X_AXIS]);
    // THEKERNEL->streams->printf("DEBUG: cartesian_mm[Y_AXIS] =  %f\n", cartesian_mm[Y_AXIS]);
    // THEKERNEL->streams->printf("DEBUG: actuator_to_cartesian (end)\n");
}

bool SCARASolution::set_optional(const arm_options_t &options)
{

    arm_options_t::const_iterator i;

    i = options.find('T'); // Theta arm1 length
    if (i != options.end())
    {
        arm1_length = i->second;
    }
    i = options.find('P'); // Psi arm2 length
    if (i != options.end())
    {
        arm2_length = i->second;
    }
    i = options.find('X'); // Home initial position X
    if (i != options.end())
    {
        scara_offset_x = i->second;
    }
    i = options.find('Y'); // Home initial position Y
    if (i != options.end())
    {
        scara_offset_y = i->second;
    }
    i = options.find('A'); // Scaling X_AXIS
    if (i != options.end())
    {
        scara_scaling_x = i->second;
    }
    i = options.find('B'); // Scaling Y_AXIS
    if (i != options.end())
    {
        scara_scaling_y = i->second;
    }
    //i= options.find('C');          // Scaling Z_AXIS
    //if(i != options.end()) {
    //    scara_scaling_z= i->second;
    //}

    init();
    return true;
}

bool SCARASolution::get_optional(arm_options_t &options, bool force_all) const
{
    options['T'] = this->arm1_length;
    options['P'] = this->arm2_length;
    options['X'] = this->scara_offset_x;
    options['Y'] = this->scara_offset_y;
    options['A'] = this->scara_scaling_x;
    options['B'] = this->scara_scaling_y;
    // options['C']= this->scara_scaling_z;

    return true;
};
