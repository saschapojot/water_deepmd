//
// Created by adada on 12/3/2025.
//
#include "./mc_read_load_compute.hpp"


void mc_computation::load_pickle_data(const std::string& filename, std::shared_ptr<double[]>& data_ptr,
                                      std::size_t size)
{
    // Initialize Python and NumPy
    Py_Initialize();
    np::initialize();


    try
    {
        // Use Python's 'io' module to open the file directly in binary mode
        py::object io_module = py::import("io");
        py::object file = io_module.attr("open")(filename, "rb"); // Open file in binary mode

        // Import the 'pickle' module
        py::object pickle_module = py::import("pickle");

        // Use pickle.load to deserialize from the Python file object
        py::object loaded_data = pickle_module.attr("load")(file);

        // Close the file
        file.attr("close")();

        // Check if the loaded object is a NumPy array
        if (py::extract<np::ndarray>(loaded_data).check())
        {
            np::ndarray np_array = py::extract<np::ndarray>(loaded_data);

            // Convert the NumPy array to a Python list using tolist()
            py::object py_list = np_array.attr("tolist")();

            // Ensure the list size matches the expected size
            ssize_t list_size = py::len(py_list);
            if (static_cast<std::size_t>(list_size) > size)
            {
                throw std::runtime_error("The provided shared_ptr array size is smaller than the list size.");
            }

            // Copy the data from the Python list to the shared_ptr array
            for (ssize_t i = 0; i < list_size; ++i)
            {
                data_ptr[i] = py::extract<double>(py_list[i]);
            }
        }
        else
        {
            throw std::runtime_error("Loaded data is not a NumPy array.");
        }
    }
    catch (py::error_already_set&)
    {
        PyErr_Print();
        throw std::runtime_error("Python error occurred.");
    }
}

void mc_computation::save_array_to_pickle(const std::shared_ptr<double[]>& ptr, int size, const std::string& filename)
{
    using namespace boost::python;
    namespace np = boost::python::numpy;

    // Initialize Python interpreter if not already initialized
    if (!Py_IsInitialized())
    {
        Py_Initialize();
        if (!Py_IsInitialized())
        {
            throw std::runtime_error("Failed to initialize Python interpreter");
        }
        np::initialize(); // Initialize NumPy
    }

    try
    {
        // Import the pickle module
        object pickle = import("pickle");
        object pickle_dumps = pickle.attr("dumps");

        // Convert C++ array to NumPy array using shared_ptr
        np::ndarray numpy_array = np::from_data(
            ptr.get(), // Use shared_ptr's raw pointer
            np::dtype::get_builtin<double>(), // NumPy data type (double)
            boost::python::make_tuple(size), // Shape of the array (1D array)
            boost::python::make_tuple(sizeof(double)), // Strides
            object() // Optional base object
        );

        // Serialize the NumPy array using pickle.dumps
        object serialized_array = pickle_dumps(numpy_array);

        // Extract the serialized data as a string
        std::string serialized_str = extract<std::string>(serialized_array);

        // Write the serialized data to a file
        std::ofstream file(filename, std::ios::binary);
        if (!file)
        {
            throw std::runtime_error("Failed to open file for writing");
        }
        file.write(serialized_str.data(), serialized_str.size());
        file.close();

        // Debug output (optional)
        // std::cout << "Array serialized and written to file successfully." << std::endl;
    }
    catch (const error_already_set&)
    {
        PyErr_Print();
        std::cerr << "Boost.Python error occurred." << std::endl;
    } catch (const std::exception& e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
}



void mc_computation::init_coord()
{

    std::string name;
    std::string coord_inFileName;
    if (this->flushLastFile == -1)
    {
        name = "init";
        coord_inFileName=out_coord_path+"/"+name+".coord.pkl";
        this->load_pickle_data(coord_inFileName,coord_init,Nx*Ny*Nz*3*3);
    }//end flushLastFile==-1
    else
    {
        name="flushEnd"+std::to_string(this->flushLastFile);
        coord_inFileName=out_coord_path+"/"+name+".coord.pkl";
        this->load_pickle_data(coord_inFileName,coord_data_all_ptr,sweep_to_write*Nx*Ny*Nz*3*3);

        //copy last Nx*Ny*Nz*3*3 elements of coord_data_all_ptr
        std::memcpy(coord_init.get(),coord_data_all_ptr.get()+(sweep_to_write-1)*Nx*Ny*Nz*3*3,Nx*Ny*Nz*3*3*sizeof(double));

    }//end else
    // print_shared_ptr(coord_init,Nx*Ny*Nz*3*3);

}


void mc_computation::init_and_run()
{

    this->init_coord();
    this->execute_mc(coord_init,newFlushNum);
}


///
/// @param x
/// @param leftEnd
/// @param rightEnd
/// @param eps
/// @return return a value within distance eps from x, on the open interval (leftEnd, rightEnd)
double mc_computation::generate_uni_open_interval(const double& x, const double& leftEnd, const double& rightEnd,
                                                  const double& eps)
{
    double xMinusEps = x - eps;
    double xPlusEps = x + eps;

    double unif_left_end = xMinusEps < leftEnd ? leftEnd : xMinusEps;
    double unif_right_end = xPlusEps > rightEnd ? rightEnd : xPlusEps;

    //    std::random_device rd;
    //    std::ranlux24_base e2(rd());

    double unif_left_end_double_on_the_right = std::nextafter(unif_left_end, std::numeric_limits<double>::infinity());


    std::uniform_real_distribution<> distUnif(unif_left_end_double_on_the_right, unif_right_end);
    //(unif_left_end_double_on_the_right, unif_right_end)

    double xNext = distUnif(e2);
    return xNext;
}


///
/// @param x proposed value
/// @param y current value
/// @param a left end of interval
/// @param b right end of interval
/// @param epsilon half length
/// @return proposal probability S(x|y)
double mc_computation::S_uni(const double& x, const double& y, const double& a, const double& b, const double& epsilon)
{
    if (a < y and y < a + epsilon)
    {
        return 1.0 / (y - a + epsilon);
    }
    else if (a + epsilon <= y and y < b + epsilon)
    {
        return 1.0 / (2.0 * epsilon);
    }
    else if (b - epsilon <= y and y < b)
    {
        return 1.0 / (b - y + epsilon);
    }
    else
    {
        std::cerr << "value out of range." << std::endl;
        std::exit(10);
    }
}



void mc_computation::coord_update( const std::vector<double>& coord_1_frame_curr,
        const std::vector<double>& coord_1_frame_next,
        double& UCurr, double& UNext)
{
//compute UCurr
    this->dp_model.compute(UCurr,force_1_frame_curr,virial_1_frame_curr,coord_1_frame_curr,atom_type,cell);

    //compute UNext
    this->dp_model.compute(UNext,force_1_frame_next,virial_1_frame_next,coord_1_frame_next,atom_type,cell);


}



double mc_computation::acceptanceRatio_uni(const std::vector<double>&coord_1_frame_curr,
        const std::vector<double>& coord_1_frame_next,
        const int& ind, const double& UCurr, const double& UNext)
{

    double numerator = -this->beta * UNext;
    double denominator = -this->beta * UCurr;
    double R = std::exp(numerator - denominator);
    double S_curr_next,S_next_curr;

    //if it is x component
    if (ind%3==0)
    {
         S_curr_next = S_uni(coord_1_frame_curr[ind],coord_1_frame_next[ind],
            0.0,box_x,h);

         S_next_curr = S_uni(coord_1_frame_next[ind],coord_1_frame_curr[ind],
            0.0,box_x,h);
    }//end mod ==0

    else if (ind%3==1)
    {
         S_curr_next = S_uni(coord_1_frame_curr[ind],coord_1_frame_next[ind],
            0.0,box_y,h);
         S_next_curr = S_uni(coord_1_frame_next[ind],coord_1_frame_curr[ind],
            0.0,box_y,h);
    }//end mod ==1

    else
    {
         S_curr_next = S_uni(coord_1_frame_curr[ind],coord_1_frame_next[ind],
            0.0,box_z,h);
         S_next_curr = S_uni(coord_1_frame_next[ind],coord_1_frame_curr[ind],
            0.0,box_z,h);
    }//end mod==2
    double ratio = S_curr_next / S_next_curr;

    if (std::fetestexcept(FE_DIVBYZERO))
    {
        std::cout << "Division by zero exception caught." << std::endl;
        std::exit(15);
    }
    if (std::isnan(ratio))
    {
        std::cout << "The result is NaN." << std::endl;
        std::exit(15);
    }

    R *= ratio;

    return std::min(1.0, R);

}


void mc_computation::execute_mc_one_sweep(std::vector<double>& coord_1_frame_curr_vec, double& UCurr,
    std::vector<double>& coord_1_frame_next_vec)
{

    double UNext = 0;

    //update coord
    for (int i=0;i<Nx*Ny*Nz*3*3;i++)
    {
        int ind=unif_in_0_NxNyNz_3_3(e2);
        // std::cout<<"ind="<<ind<<std::endl;
        //proposal
        this->proposal_uni(coord_1_frame_curr_vec,coord_1_frame_next_vec,ind);
        //energy
        this->coord_update(coord_1_frame_curr_vec,coord_1_frame_next_vec,UCurr,UNext);
        //acceptance ratio
        double r = this->acceptanceRatio_uni(coord_1_frame_curr_vec,coord_1_frame_next_vec,ind,
            UCurr,UNext);
        double u = distUnif01(e2);
        // std::cout<<"u="<<u<<std::endl;
        if (u<=r)
        {
            UCurr=UNext;
            std::memcpy(coord_1_frame_curr_vec.data(),coord_1_frame_next_vec.data(),Nx*Ny*Nz*3*3*sizeof(double));
        }//end of accept-reject

    }//end updating coord


}


void mc_computation::proposal_uni(const std::vector<double>&coord_1_frame_curr,
       std::vector<double>& coord_1_frame_next,const int& ind )
{
    // size_t length=coord_1_frame_curr.size();

    double elem_val_new=0;
    if (ind%3==0)
    {
        elem_val_new=this->generate_uni_open_interval(coord_1_frame_curr[ind],
            0.0,box_x,h);
    }//end mod ==0
    else if (ind%3==1)
    {
        elem_val_new=this->generate_uni_open_interval(coord_1_frame_curr[ind],
            0.0,box_y,h);
    }//end mod ==1
    else
    {
        elem_val_new=this->generate_uni_open_interval(coord_1_frame_curr[ind],
            0.0,box_z,h);
    }//end mod ==2

    std::memcpy(coord_1_frame_next.data(),coord_1_frame_curr.data(),Nx*Ny*Nz*3*3*sizeof(double));
    coord_1_frame_next[ind]=elem_val_new;
}


void mc_computation::execute_mc(const std::shared_ptr<double[]>& coord_1_frame,const int& flushNum)
{
std::memcpy(coord_1_frame_curr.data(),coord_1_frame.get(),Nx*Ny*Nz*3*3*sizeof(double));

    // print_vec(coord_1_frame_curr,Nx*Ny*Nz*3*3);
    // print_vec(coord_1_frame_next,Nx*Ny*Nz*3*3);

    double UCurr=0;
    int flushThisFileStart=this->flushLastFile+1;
    for (int fls=0;fls<flushNum;fls++)
    {
        const auto tMCStart{std::chrono::steady_clock::now()};
        for (int swp = 0; swp < sweep_to_write*sweep_multiple; swp++)
        {

            this->execute_mc_one_sweep(coord_1_frame_curr,UCurr,coord_1_frame_next);


            if (swp%sweep_multiple==0)
            {
                int swp_out=swp/sweep_multiple;

                this->U_data_all_ptr[swp_out]=UCurr;
                std::memcpy(coord_data_all_ptr.get()+swp_out*Nx*Ny*Nz*3*3,
                    coord_1_frame_curr.data(),Nx*Ny*Nz*3*3*sizeof(double));
            }//end save to array

        }//end sweep for

        int flushEnd=flushThisFileStart+fls;
        std::string fileNameMiddle =  "flushEnd" + std::to_string(flushEnd);

        std::string out_U_PickleFileName = out_U_path+"/" + fileNameMiddle + ".U.pkl";

        std::string out_coord_PickleFileName=out_coord_path+"/"+fileNameMiddle+".coord.pkl";

        //save U
        this->save_array_to_pickle(U_data_all_ptr,sweep_to_write,out_U_PickleFileName);

        //save coord
        this->save_array_to_pickle(coord_data_all_ptr,sweep_to_write*Nx*Ny*Nz*3*3,out_coord_PickleFileName);

        const auto tMCEnd{std::chrono::steady_clock::now()};
        const std::chrono::duration<double> elapsed_secondsAll{tMCEnd - tMCStart};
        std::cout << "flush " + std::to_string(flushEnd)  + ": "
                  << elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;
    }//end fls
    std::cout << "mc executed for " << flushNum << " flushes." << std::endl;

}
