//
// Created by adada on 12/3/2025.
//

#ifndef MC_READ_LOAD_COMPUTE_HPP
#define MC_READ_LOAD_COMPUTE_HPP
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <cfenv> // for floating-point exceptions
#include "deepmd/DeepPot.h"
#include <fstream>
#include <Python.h>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#endif //MC_READ_LOAD_COMPUTE_HPP

namespace fs = boost::filesystem;
namespace py = boost::python;
namespace np = boost::python::numpy;
static const int seed=29;
class mc_computation
{
    public:
    mc_computation(const std::string& cppInParamsFileName): e2(seed), distUnif01(0.0, 1.0)
    {

        std::ifstream file(cppInParamsFileName);
        if (!file.is_open())
        {
            std::cerr << "Failed to open the file." << std::endl;
            std::exit(20);
        }
        std::string line;
        int paramCounter = 0;
        while (std::getline(file, line))
        {
            // Check if the line is empty
            if (line.empty())
            {
                continue; // Skip empty lines
            }
            std::istringstream iss(line);

            //read T
            if (paramCounter == 0)
            {
                iss >> T;
                if (T <= 0)
                {
                    std::cerr << "T must be >0" << std::endl;
                    std::exit(1);
                } //end if
                std::cout << "T=" << T <<"K"<< std::endl;
                this->beta = 1.0 / T*1.0/kB;
                std::cout << "beta=" << beta << std::endl;
                paramCounter++;
                continue;
            } //end T
            //read P
            if (paramCounter == 1)
            {
                iss>>pressure;
                std::cout<<"pressure="<<pressure<<"GPa"<<std::endl;
                this->gamma=1.0/1.380649*1e2*pressure/T;
                std::cout<<"gamma="<<gamma<<std::endl;
                paramCounter++;
                continue;
            }

            //read model_file
            if (paramCounter == 2)
            {
                iss>>model_file;
                dp_model=deepmd::DeepPot(model_file);
                paramCounter++;
                continue;
            }//end  model_file

            //read sweep_to_write
            if (paramCounter == 3)
            {
                iss>>sweep_to_write;
                if (sweep_to_write<=0)
                {
                    std::cerr << "sweep_to_write must be >0" << std::endl;
                    std::exit(1);
                }
                std::cout << "sweep_to_write=" << sweep_to_write << std::endl;
                paramCounter++;
                continue;
            }//end sweep_to_write

            //read flushLastFile
            if (paramCounter == 4)
            {
                iss>>flushLastFile;
                std::cout << "flushLastFile=" << flushLastFile << std::endl;
                paramCounter++;
                continue;

            }//end flushLastFile

            //newFlushNum
            if (paramCounter == 5)
            {
                iss>>newFlushNum;
                if (newFlushNum <= 0)
                {
                    std::cerr << "newFlushNum must be >0" << std::endl;
                    std::exit(1);
                }
                std::cout << "newFlushNum=" << newFlushNum << std::endl;
                paramCounter++;
                continue;
            }//end newFlushNum

        // read TDirRoot
            if (paramCounter == 6)
            {
                iss>>TDirRoot;
                std::cout<<"TDirRoot="<<TDirRoot<<std::endl;
                paramCounter++;
                continue;
            }//end TDirRoot

            //read U_dist_dataDir
            if (paramCounter == 7)
            {
                iss>>U_dist_dataDir;
                std::cout<<"U_dist_dataDir="<<U_dist_dataDir<<std::endl;
                paramCounter++;
                continue;
            }//end U_dist_dataDir

            //read h
            if (paramCounter == 8)
            {
                iss>>h;
                if (h <= 0)
                {
                    std::cerr << "h must be >0" << std::endl;
                    std::exit(1);
                }
                std::cout << "h=" << h << std::endl;
                paramCounter++;
                continue;
            }//end h

            //read sweep_multiple
            if (paramCounter == 9)
            {
                iss>>sweep_multiple;
                if (sweep_multiple <= 0)
                {
                    std::cerr << "sweep_multiple must be >0" << std::endl;
                    std::exit(1);
                }
                std::cout << "sweep_multiple=" << sweep_multiple << std::endl;
                paramCounter++;
                continue;
            }//end sweep_multiple

            //read Nx
            if (paramCounter == 10)
            {
                iss>>Nx;
                if (Nx<=0)
                {
                    std::cerr << "Nx must be >0" << std::endl;
                    std::exit(1);
                }
                std::cout<<"Nx="<<Nx<<std::endl;
                paramCounter++;
                continue;
            }//end Nx

            //read Ny
            if (paramCounter == 11)
            {
                iss>>Ny;
                if (Ny<=0)
                {
                    std::cerr << "Ny must be >0" << std::endl;
                    std::exit(1);
                }
                std::cout<<"Ny="<<Ny<<std::endl;
                paramCounter++;
                continue;

            }//end Ny

            //read Nz
            if (paramCounter == 12)
            {
                iss>>Nz;
                if (Nz<=0)
                {
                    std::cerr << "Nz must be >0" << std::endl;
                    std::exit(1);
                }
                std::cout<<"Nz="<<Nz<<std::endl;
                paramCounter++;
                continue;

            }//end  Nz

            // read box_upper_bound
            if (paramCounter == 13)
            {
                iss>>box_upper_bound;
                std::cout<<"box_upper_bound="<<box_upper_bound<<std::endl;

                paramCounter++;
                continue;
            }//end box_upper_bound






        }//end getline while


        this->total_atom_num=Nx*Ny*Nz*3;
        this->total_atom_xyz_components_num=Nx*Ny*Nz*3*3;
        this->box_component_num=9;
        this->box_x_component_position=0;
        this->box_y_component_position=4;
        this->box_z_component_position=8;
        non_0_inds_in_box={box_x_component_position,box_y_component_position,box_z_component_position};
        //allocate memory for data

        try
        {
            this->U_data_all_ptr = std::shared_ptr<double[]>(new double[sweep_to_write],
                                                             std::default_delete<double[]>());
            this->coord_data_all_ptr=std::shared_ptr<double[]>(new double[sweep_to_write*total_atom_xyz_components_num],
                                                             std::default_delete<double[]>());

            this->box_data_all_ptr=std::shared_ptr<double[]>(new double[sweep_to_write*box_component_num],
                                                             std::default_delete<double[]>());
            this->coord_init=std::shared_ptr<double[]>(new double[total_atom_xyz_components_num],
                                                             std::default_delete<double[]>());
            this->coord_init_vector.resize(total_atom_xyz_components_num);

            this->box_init=std::shared_ptr<double[]>(new double[box_component_num],
                                                             std::default_delete<double[]>());

            this->box_init_vector.resize(box_component_num);

            this->coord_x_components.resize(total_atom_num);

            this->coord_y_components.resize(total_atom_num);

            this->coord_z_components.resize(total_atom_num);


            this->coord_1_frame_curr.resize(total_atom_xyz_components_num);
            this->coord_1_frame_next.resize(total_atom_xyz_components_num);

            this->box_1_frame_curr.resize(box_component_num);
            this->box_1_frame_next.resize(box_component_num);

            this->force_1_frame_curr.resize(total_atom_xyz_components_num);
            this->force_1_frame_next.resize(total_atom_xyz_components_num);

            this->virial_1_frame_curr.resize(total_atom_xyz_components_num);
            this->virial_1_frame_next.resize(total_atom_xyz_components_num);


            this->atom_type.resize(total_atom_num);

        }
        catch (const std::bad_alloc& e)
        {
            std::cerr << "Memory allocation error: " << e.what() << std::endl;
            std::exit(2);
        } catch (const std::exception& e)
        {
            std::cerr << "Exception: " << e.what() << std::endl;
            std::exit(2);
        }



        // print_vec(cell,9);

        // print_shared_ptr(box_data_ptr,9);
        this->type_file=U_dist_dataDir+"/raw.pkl";

        std::shared_ptr<double[]> type_data_ptr=std::shared_ptr<double[]>(new double[total_atom_num],
                                                             std::default_delete<double[]>());

        this->load_pickle_data(type_file,type_data_ptr,total_atom_num);
        // print_shared_ptr(type_data_ptr,10);
        // std::cout<<type_data_ptr[1]<<std::endl;
        for (int i=0;i<total_atom_num;i++)
        {
            this->atom_type[i]=static_cast<int>( type_data_ptr[i]);
        }//end for
        // print_vec(atom_type,Nx*Ny*Nz*3);

        this->unif_in_0_NxNyNz_3_3_m1=std::uniform_int_distribution<int>(0,total_atom_xyz_components_num-1);

        this->unif_in_0_2=std::uniform_int_distribution<int>(0,2);

        // this->box_upper_bound=15.0;
        //


        this->out_U_path=U_dist_dataDir+"/U/";
        this->out_coord_path=U_dist_dataDir+"/coord/";
        this->out_box_path=U_dist_dataDir+"/box/";
    }//end constructor

public:
    void init_and_run();

    void init_coord_and_box();

    void execute_mc(const std::shared_ptr<double[]>&coord_init,
        const std::shared_ptr<double[]>&box_init,const int& flushNum);

    void execute_mc_one_sweep(std::vector<double>& coord_1_frame_curr,
        std::vector<double>&box_1_frame_curr,
        double& UCurr,
       std::vector<double>&coord_1_frame_next,
       std::vector<double>&box_1_frame_next);


    void energy_update_coord_one_component(const std::vector<double>& coord_1_frame_curr,
        const std::vector<double>&coord_1_frame_next,  double& UCurr,   double& UNext,const std::vector<double>&box_1_frame_curr);


    void energy_update_box_one_component(const std::vector<double>&box_1_frame_curr,
        const std::vector<double>&box_1_frame_next,const std::vector<double>&coord_1_frame_curr,double& UCurr, double& UNext);

    double acceptanceRatio_uni_for_coord(const std::vector<double>& coord_1_frame_curr,
        const std::vector<double>&coord_1_frame_next,const int& ind,const double& UCurr, const double& UNext,
       const double& box_x_val, const double& box_y_val,const double& box_z_val );
    ///
    /// @param box_1_frame_curr
    /// @param box_1_frame_next
    /// @param ind
    /// @param UCurr
    /// @param UNext
    /// @return acceptance ratio for updating box
    double acceptanceRatio_uni_for_box(const std::vector<double>&box_1_frame_curr,
                                         const std::vector<double>&box_1_frame_next,const int& ind,
                                         const double& UCurr, const double& UNext,
                                         const double &atom_x_max_val, const double &atom_y_max_val, const double &atom_z_max_val);


    void proposal_uni_box(const std::vector<double>&box_1_frame_curr,
        std::vector<double>&box_1_frame_next,const int& ind,const std::vector<double>&coord_1_frame_curr,
        double &atom_x_max_val,double & atom_y_max_val,double & atom_z_max_val);

    void proposal_uni_coord(const std::vector<double>& coord_1_frame_curr,
        std::vector<double>&coord_1_frame_next,const int& ind,const std::vector<double>&box_1_frame_curr);


    double box_2_volume(const std::vector<double>& box_1_frame);
    ///
    /// @param x proposed value
    /// @param y current value
    /// @param a left end of interval
    /// @param b right end of interval
    /// @param epsilon half length
    /// @return proposal probability S(x|y)
    double S_uni(const double& x, const double& y, const double& a, const double& b, const double& epsilon);
    ///
    /// @param x
    /// @param leftEnd
    /// @param rightEnd
    /// @param eps
    /// @return return a value within distance eps from x, on the open interval (leftEnd, rightEnd)
    double generate_uni_open_interval(const double& x, const double& leftEnd, const double& rightEnd,
                                      const double& eps);

    void load_pickle_data(const std::string& filename, std::shared_ptr<double[]>& data_ptr, std::size_t size);
    void save_array_to_pickle(const std::shared_ptr<double[]>& ptr, int size, const std::string& filename);

    template <class T>
    void print_shared_ptr(const std::shared_ptr<T>& ptr, const int& size)
    {
        if (!ptr)
        {
            std::cout << "Pointer is null." << std::endl;
            return;
        }

        for (int i = 0; i < size; i++)
        {
            if (i < size - 1)
            {
                std::cout << ptr[i] << ",";
            }
            else
            {
                std::cout << ptr[i] << std::endl;
            }
        }
    } //end print_shared_ptr

    template <class T>
    void print_vec(const std::vector<T>& vec, const int& size)
    {
        if (vec.empty())
        {
            std::cout << "Vector is empty." << std::endl;
            return;
        }
        for (int i = 0; i < size; i++)
        {
            if (i < size - 1)
            {
                std::cout << vec[i] << ",";
            }
            else
            {
                std::cout << vec[i] << std::endl;
            }
        }//end for
    }//end print_vec

public:
    double T; // temperature
    double beta;
    double pressure;
    double gamma;
    double h; // step size
    int sweep_to_write;
    int newFlushNum;
    int flushLastFile;
    std::string TDirRoot;
    std::string U_dist_dataDir;
    std::ranlux24_base e2;
    std::uniform_int_distribution<int> unif_in_0_NxNyNz_3_3_m1;
    std::uniform_int_distribution<int> unif_in_0_2;
    std::uniform_real_distribution<> distUnif01;
    int sweep_multiple;
    std::string out_U_path;
    std::string out_coord_path;
    std::string out_box_path;
    std::string model_file;

    double kB=8.617333262e-5;

    deepmd::DeepPot dp_model;
    int Nx,Ny,Nz;

    double box_upper_bound;

    //data in 1 flush
    std::shared_ptr<double[]> U_data_all_ptr; //all U data
    std::shared_ptr<double[]> coord_data_all_ptr;// all coord data, O,H,H,O,H,H,...,O,H,H
    std::shared_ptr<double[]> box_data_all_ptr;//all box data
    //initial value
    std::shared_ptr<double[]> coord_init;
    std::shared_ptr<double[]>box_init;

    std::vector<double>coord_init_vector;
    std::vector<double>box_init_vector;

    std::vector<double>coord_1_frame_curr;
    std::vector<double>coord_1_frame_next;
    std::vector<double> box_1_frame_curr;
    std::vector<double> box_1_frame_next;

    std::vector<int>atom_type;

    std::vector<double>force_1_frame_curr,virial_1_frame_curr;
    std::vector<double>force_1_frame_next,virial_1_frame_next;


    std::vector<double> coord_x_components,coord_y_components,coord_z_components;

    std::string type_file;

    //for 1 frame:
    int total_atom_num;
    int total_atom_xyz_components_num;

    int box_component_num;

    int box_x_component_position,box_y_component_position,box_z_component_position;

    std::vector<int> non_0_inds_in_box;

};