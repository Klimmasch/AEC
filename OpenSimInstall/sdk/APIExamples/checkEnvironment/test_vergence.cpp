#include <OpenSim/OpenSim.h>

using namespace OpenSim;
using namespace SimTK;

//______________________________________________________________________________
/**
 * Create a model that does nothing. 
 */
int main(int argc, char *argv[])
{
    double rotation_x = 0;
    double rotation_y = 0;
    double rotation_z = 0;
	try {
        if (argc < 6)
        {
            std::cerr << "NOT ENOUGH INPUT PARAMS" << std::endl;
            return 1;
        }
        if (argc > 6)
        {
            std::cerr << "TOO MANY INPUT PARAMS" << std::endl;
            return 1;
        }
        if (argc == 6)
        {

            std::string texture = std::string(argv[1]);
            //std::string texture_right = std::string(argv[2]);
            float distance = atof(std::string(argv[2]).c_str());
            //float distance_right = atof(std::string(argv[4]).c_str());
            //rotation_x = atof(std::string(argv[3]).c_str());
            rotation_y = atof(std::string(argv[3]).c_str());
            //rotation_z = atof(std::string(argv[5]).c_str());
            std::string save_left = std::string(argv[4]);
            std::string save_right = std::string(argv[5]);

            //rotation = atof(std::string(argv[7]).c_str()); // Here we have a vergence angle (one between 2 lines of signs)

            ///////////////////////////////////////////
            // DEFINE BODIES AND JOINTS OF THE MODEL //
            ///////////////////////////////////////////
            Model osimModel("eye_model.osim");


            Vec3 left(0);
            left[0] = -0.0307;
            left[1] = 1.578;
            left[2] = -0.028;
            Vec3 right(0);
            right[0] = -0.0307;
            right[1] = 1.578;
            right[2] = 0.028;

            //double x2Rot = rotation_x*(3.141617/180);
            double x2Rot = 0;
            double y2Rot = (rotation_y/2)*(3.141617/180);
            double z2Rot = 0;
            //double z2Rot = rotation_z*(3.141617/180);



            Rotation leftRot;
            leftRot.setRotationFromThreeAnglesThreeAxes(SimTK::BodyRotationSequence,x2Rot,XAxis,-y2Rot-1.57079633,YAxis,z2Rot,ZAxis);
            Rotation rightRot;
            rightRot.setRotationFromThreeAnglesThreeAxes(SimTK::BodyRotationSequence,x2Rot,XAxis,y2Rot-1.57079633,YAxis,z2Rot,ZAxis);

            osimModel.buildSystem();
            ModelVisualizer* viz2 = new ModelVisualizer(osimModel,true);

            SimTK::State& si = osimModel.initializeState(viz2);
            osimModel.getMultibodySystem().realize(si, Stage::Position);
            // Compute initial conditions for muscles
            osimModel.equilibrateMuscles(si);



            // start timing
            std::clock_t startTime = std::clock();

            viz2->updSimbodyVisualizer().setTexture(texture);
            viz2->updSimbodyVisualizer().initRenderer();
            viz2->updSimbodyVisualizer().setDistance(distance);
            viz2->updSimbodyVisualizer().setCameraTransformLeft(Transform(leftRot,left));
            viz2->updSimbodyVisualizer().setCameraTransformRight(Transform(rightRot,right));
            viz2->updSimbodyVisualizer().setImageLeft(save_left);
            viz2->updSimbodyVisualizer().setImageRight(save_right);

            //std::cout << "*" << rotation_y << "*";

            viz2->updSimbodyVisualizer().drawFrameNow(si,true);

            viz2->updSimbodyVisualizer().rendererScene();
//            std::cout << "Time";
            std::cout<< "Elapsed time:" << (std::clock() - startTime )/(CLOCKS_PER_SEC/1000) << "s" << "\n" ;
        }
	}
	catch (const std::exception& ex)
    {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr << "UNRECOGNIZED EXCEPTION" << std::endl;
        return 1;
    }

    std::cout << "OpenEyeSim test completed successfully. " << std::endl;

	return 0;
}
