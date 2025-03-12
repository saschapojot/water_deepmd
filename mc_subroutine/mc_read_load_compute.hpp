//
// Created by adada on 12/3/2025.
//

#ifndef MC_READ_LOAD_COMPUTE_HPP
#define MC_READ_LOAD_COMPUTE_HPP

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

class mc_computation
{
    public:
    mc_computation(const std::string& cppInParamsFileName): e2(std::random_device{}()), distUnif01(0.0, 1.0)
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

            }

            // read box_x
            if (paramCounter == 13)
            {
                iss>>box_x;
                if (box_x<=0)
                {
                    std::cerr << "box_x must be >0" << std::endl;
                    std::exit(1);
                }
                std::cout<<"box_x="<<box_x<<std::endl;
                paramCounter++;
                continue;
            }//end box_x

            //read box_y
            if (paramCounter == 14)
            {
                iss>>box_y;
                if (box_y<=0)
                {
                    std::cerr << "box_y must be >0" << std::endl;
                    std::exit(1);
                }
                std::cout<<"box_y="<<box_y<<std::endl;
                paramCounter++;
                continue;

            }//end box_y

            //read box_z
            if (paramCounter == 15)
            {
                iss>>box_z;
                if (box_z<=0)
                {
                    std::cerr << "box_z must be >0" << std::endl;
                    std::exit(1);
                }
                std::cout<<"box_z="<<box_z<<std::endl;
                paramCounter++;
                continue;
            }//end box_z
        }//end getline while

        //allocate memory for data

        try
        {
            this->U_data_all_ptr = std::shared_ptr<double[]>(new double[sweep_to_write],
                                                             std::default_delete<double[]>());
            this->coord_data_all_ptr=std::shared_ptr<double[]>(new double[sweep_to_write*Nx*Ny*Nz*3*3],
                                                             std::default_delete<double[]>());

            this->coord_1_frame.reserve(Nx*Ny*Nz*3*3);
            this->force_1_frame.reserve(Nx*Ny*Nz*3*3);

            this->virial_1_frame.reserve(Nx*Ny*Nz*3*3);

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
        this->out_U_path = this->U_dist_dataDir + "/U/";
        if (!fs::is_directory(out_U_path) || !fs::exists(out_U_path))
        {
            fs::create_directories(out_U_path);
        }

        this->out_coord_path=this->U_dist_dataDir + "/coord/";
        if (!fs::is_directory(out_coord_path) || !fs::exists(out_coord_path))
        {
            fs::create_directories(out_coord_path);
        }

        //read box
        this->box_file=U_dist_dataDir+"/box.pkl";

        std::shared_ptr<double[]> box_data_ptr=std::shared_ptr<double[]>(new double[9],
                                                             std::default_delete<double[]>());
        this->load_pickle_data(box_file,box_data_ptr,9);


        this->type_file=U_dist_dataDir+"/raw.pkl";

        std::shared_ptr<double[]> type_data_ptr=std::shared_ptr<double[]>(new double[Nx*Ny*Nz*3],
                                                             std::default_delete<double[]>());

        // this->load_pickle_data(type_file,type_data_ptr,Nx*Ny*Nz*3);
        // print_shared_ptr(type_data_ptr,10);
    }//end constructor

public:
    void init_and_run();
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
    std::uniform_real_distribution<> distUnif01;
    int sweep_multiple;
    std::string out_U_path;
    std::string out_coord_path;
    std::string model_file;

    double kB=8.617333262e-5;

    deepmd::DeepPot dp_model;
    int Nx,Ny,Nz;
    double box_x,box_y,box_z;
    //data in 1 flush
    std::shared_ptr<double[]> U_data_all_ptr; //all U data
    std::shared_ptr<double[]> coord_data_all_ptr;// all coord data, O,H,H,O,H,H,...,O,H,H

    std::vector<double>coord_1_frame;
    std::vector<int>atom_type;
    std::vector<double> cell;
    std::vector<double>force_1_frame,virial_1_frame;

    std::string box_file;
    std::string type_file;

};