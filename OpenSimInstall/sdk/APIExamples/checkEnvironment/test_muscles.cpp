
#include <OpenSim/OpenSim.h>

using namespace OpenSim;
using namespace SimTK;

//______________________________________________________________________________
/**
 * Create a model that does nothing.
 */
int main(int argc, char *argv[])
{
    try {
        if (argc < 5)
        {
            std::cerr << "NOT ENOUGH INPUT PARAMS" << std::endl;
            return 1;
        }
        if (argc > 5)
        {
            std::cerr << "TOO MANY INPUT PARAMS" << std::endl;
            return 1;
        }
        if (argc == 5)
        {
            double rotation = 0;

            std::string texture = std::string(argv[1]);
            float distance = atof(std::string(argv[2]).c_str());
            std::string save_left = std::string(argv[3]);
            std::string save_right = std::string(argv[4]);


            std::cout<<"Constructing tool from setup file "<<".\n\n";
            ForwardTool forward("eye_Setup_Forward.xml");
            // PRINT MODEL INFORMATION
            Model& osimModel = forward.getModel();
            Manager manager = forward.runManager();
            Array<SimTK::State> states = manager.savedStates;
            Storage simulated = manager.getStateStorage();

            Array <double> x2Rot,y2Rot,z2Rot,x2Tra,y2Tra,z2Tra;
            simulated.getDataColumn("left_xRotation",x2Rot);
            simulated.getDataColumn("left_yRotation",y2Rot);
            simulated.getDataColumn("left_zRotation",z2Rot);
            simulated.getDataColumn("left_xTranslation",x2Tra);
            simulated.getDataColumn("left_yTranslation",y2Tra);
            simulated.getDataColumn("left_zTranslation",z2Tra);

            Vec3 left(0);
            left[0] = x2Tra[0]-0.0307;
            left[1] = x2Tra[0]+1.578;
            left[2] = x2Tra[0]-0.028;
            Vec3 right(0);
            right[0] = x2Tra[0]-0.0307;
            right[1] = y2Tra[0]+1.578;
            right[2] = z2Tra[0]+0.028;

            rotation = y2Rot[states.getSize()-1];

            Rotation leftRot;
            leftRot.setRotationFromThreeAnglesThreeAxes(SimTK::BodyRotationSequence,x2Rot[states.getSize()-1],XAxis,y2Rot[states.getSize()-1]-1.57079633,YAxis,z2Rot[states.getSize()-1],ZAxis);
            Rotation rightRot;
            rightRot.setRotationFromThreeAnglesThreeAxes(SimTK::BodyRotationSequence,x2Rot[states.getSize()-1],XAxis,-y2Rot[states.getSize()-1]-1.57079633,YAxis,z2Rot[states.getSize()-1],ZAxis);


            osimModel.buildSystem();
            ModelVisualizer* viz2 = new ModelVisualizer(osimModel,true);
            osimModel.initializeState(viz2);

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
            viz2->updSimbodyVisualizer().drawFrameNow(states[states.getSize()-1],true);
            //viz2->updSimbodyVisualizer().drawFrameNow(si,true);

            viz2->updSimbodyVisualizer().rendererScene();
            std::cout << "Time";
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
