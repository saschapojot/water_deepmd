
# include "mc_read_load_compute.hpp"

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



void mc_computation::init_coord_and_box()
{
    std::string name;
    std::string coord_inFileName, box_inFileName;
    if (this->flushLastFile == -1)
    {

        name = "init";
        coord_inFileName=out_coord_path+"/"+name+".coord.pkl";

        box_inFileName=out_box_path+"/"+name+".box.pkl";

        //load coord
        this->load_pickle_data(coord_inFileName,coord_init,total_atom_xyz_components_num);

        //load box
        this->load_pickle_data(box_inFileName,box_init,box_component_num);

    }//end flushLastFile==-1
    else
    {

        name="flushEnd"+std::to_string(this->flushLastFile);

        coord_inFileName=out_coord_path+"/"+name+".coord.pkl";
        box_inFileName=out_box_path+"/"+name+".box.pkl";

        //load coord
        this->load_pickle_data(coord_inFileName,coord_data_all_ptr,
            sweep_to_write*total_atom_xyz_components_num);

        std::memcpy(coord_init.get(),
            coord_data_all_ptr.get()+(sweep_to_write-1)*total_atom_xyz_components_num,
            total_atom_xyz_components_num*sizeof(double));

        //load box
        this->load_pickle_data(box_inFileName,box_data_all_ptr,sweep_to_write*box_component_num);

        std::memcpy(box_init.get(),
            box_data_all_ptr.get()+(sweep_to_write-1)*box_component_num,box_component_num*sizeof(double));
    }//end else

    //smart pointer to std::vector, coord
    std::memcpy(coord_init_vector.data(),coord_init.get(),total_atom_xyz_components_num*sizeof(double));

    //smart pointer to std::vector, box
    std::memcpy(box_init_vector.data(),box_init.get(),box_component_num*sizeof(double));

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
    else if (a + epsilon <= y and y < b - epsilon)
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


void mc_computation::init_and_run()
{
    this->init_coord_and_box();
    // print_shared_ptr(this->coord_init,Nx*Ny*Nz*3*3);
// print_vec(this->coord_init_vector,Nx*Ny*Nz*3*3);
     // print_vec(this->box_init_vector,9);
}


void mc_computation::proposal_uni_coord(const std::vector<double>& coord_1_frame_curr,
        std::vector<double>&coord_1_frame_next,const int& ind,const std::vector<double>&box_1_frame_curr)
{
    std::memcpy(coord_1_frame_next.data(),coord_1_frame_curr.data(),total_atom_xyz_components_num*sizeof(double));

    double box_x=box_1_frame_curr[box_x_component_position];
    double box_y=box_1_frame_curr[box_y_component_position];
    double box_z=box_1_frame_curr[box_z_component_position];

    if (ind%3==0)
    {
        double x_tmp= this->generate_uni_open_interval(coord_1_frame_curr[ind],0,box_x,h);
        coord_1_frame_next[ind]=x_tmp;
    }// end mod 3 ==0
    else if (ind%3==1)
    {
        double y_tmp=this->generate_uni_open_interval(coord_1_frame_curr[ind],0,box_y,h);
        coord_1_frame_next[ind]=y_tmp;
    } //end mod 3 ==1
    else
    {
        double z_tmp=this->generate_uni_open_interval(coord_1_frame_curr[ind],0,box_z,h);
        coord_1_frame_next[ind]=z_tmp;
    }//end mod 3==2

}


void mc_computation::proposal_uni_box(const std::vector<double>&box_1_frame_curr,
        std::vector<double>&box_1_frame_next,const int& ind,const std::vector<double>&coord_1_frame_curr,
        double &atom_x_max_val,double & atom_y_max_val,double & atom_z_max_val)
{
    std::memcpy(box_1_frame_next.data(),box_1_frame_curr.data(),box_component_num*sizeof(double));
    ///////////////////////////////////////////////////////
    /// x direction of box
    if (ind==box_x_component_position)
    {
        //update box length in x direction
        //extract x components
        for (int i=0;i<total_atom_xyz_components_num;i+=3)
        {
            coord_x_components[i/3]=coord_1_frame_curr[i];
        }//end for
        auto x_max_iter=std::max_element(coord_x_components.begin(),coord_x_components.end());
        atom_x_max_val=*x_max_iter;
        //to prevent unwanted usage
        atom_y_max_val=-1;
        atom_z_max_val=-1;
        double x_direction_new_val=this->generate_uni_open_interval(box_1_frame_curr[box_x_component_position],atom_x_max_val,box_upper_bound,h);
        box_1_frame_next[box_x_component_position]=x_direction_new_val;
        return;
    }//end if ind==box_x_component_position
    /// end x direction of box
    ///////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////
    /// y direction of box
    if (ind==box_y_component_position)
    {
        //update box length in y direction
        //extract y components
        for (int i=1;i<total_atom_xyz_components_num;i+=3)
        {
            coord_y_components[(i-1)/3]=coord_1_frame_curr[i];
        }//end for
        auto max_y_iter=std::max_element(coord_y_components.begin(),coord_y_components.end());
        atom_y_max_val=*max_y_iter;
        //to prevent unwanted usage
        atom_x_max_val=-1;
        atom_z_max_val=-1;
        double y_direction_new_val=this->generate_uni_open_interval(box_1_frame_curr[box_y_component_position],atom_y_max_val,box_upper_bound,h);
        box_1_frame_next[box_y_component_position]=y_direction_new_val;
        return;
    }//end if ind==box_y_component_position
    /// end y direction of box
    ///////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////
    /// z direction of box
    if (ind==box_z_component_position)
    {
        //update box length in z direction
        //extract z components
        for (int i=2;i<total_atom_xyz_components_num;i+=3)
        {
            coord_z_components[(i-2)/3]=coord_1_frame_curr[i];
        }//end for
        auto max_z_iter=std::max_element(coord_z_components.begin(),coord_z_components.end());
        atom_z_max_val=*max_z_iter;
        //to prevent unwanted usage
        atom_x_max_val=-1;
        atom_y_max_val=-1;
        double z_direction_new_val=this->generate_uni_open_interval(box_1_frame_curr[box_z_component_position],atom_z_max_val,box_upper_bound,h);
        box_1_frame_next[box_z_component_position]=z_direction_new_val;
        return;
    }//end if ind==box_z_component_position
    /// end z direction of box
    ///////////////////////////////////////////////////////
    std::cerr<<"Invalid ind value for box: ind="<<ind<<std::endl;
    std::exit(3);
}



///
/// @param box_1_frame_curr
/// @param box_1_frame_next
/// @param ind
/// @param UCurr
/// @param UNext
/// @return acceptance ratio for updating box
double mc_computation::acceptanceRatio_uni_for_box(const std::vector<double>&box_1_frame_curr,
                                         const std::vector<double>&box_1_frame_next,const int& ind,
                                         const double& UCurr, const double& UNext,
                                         const double &atom_x_max_val, const double &atom_y_max_val, const double &atom_z_max_val)
{

    double V_next=this->box_2_volume(box_1_frame_next);
    double V_curr=this->box_2_volume(box_1_frame_curr);

    double numerator = -this->beta * UNext-this->gamma*V_next;
    double denominator = -this->beta * UCurr-this->gamma*V_curr ;
    double R = std::exp(numerator - denominator);
    double S_curr_next=0;
    double S_next_curr=0;


    if (ind==box_x_component_position)
    {
        S_curr_next=S_uni(box_1_frame_curr[ind],
            box_1_frame_next[ind],atom_x_max_val,box_upper_bound,h);

        S_next_curr=S_uni(box_1_frame_next[ind],box_1_frame_curr[ind],
           atom_x_max_val,box_upper_bound,h );
    }//end if ind==box_x_component_position
    else if (ind==box_y_component_position)
    {
        S_curr_next=S_uni(box_1_frame_curr[ind],box_1_frame_next[ind],
            atom_y_max_val,box_upper_bound,h);

        S_next_curr=S_uni(box_1_frame_next[ind],box_1_frame_curr[ind],
            atom_y_max_val,box_upper_bound,h);
    }//end else if ind==box_y_component_position

    else if (ind==box_z_component_position)
    {
        S_curr_next=S_uni(box_1_frame_curr[ind],box_1_frame_next[ind],
        atom_z_max_val,box_upper_bound,h
            );
        S_next_curr=S_uni(box_1_frame_next[ind],box_1_frame_curr[ind],
            atom_z_max_val,box_upper_bound,h);
    }//end else if ind==box_z_component_position
    else
    {
        std::cerr<<"Invalid ind="<<ind<<std::endl;
        std::exit(4);
    }//end else
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



double mc_computation::acceptanceRatio_uni_for_coord(const std::vector<double>& coord_1_frame_curr,
        const std::vector<double>&coord_1_frame_next,const int& ind,const double& UCurr, const double& UNext,
       const double& box_x_val, const double& box_y_val,const double& box_z_val )
{

    double numerator = -this->beta * UNext;
    double denominator = -this->beta * UCurr;
    double R = std::exp(numerator - denominator);
    double S_curr_next=0;
    double S_next_curr=0;


    if (ind%3==0)
    {
        S_curr_next=S_uni(coord_1_frame_curr[ind],coord_1_frame_next[ind],
            0,box_x_val,h);
        S_next_curr=S_uni(coord_1_frame_next[ind],coord_1_frame_curr[ind],
            0,box_x_val,h);
    }//end mod 3 ==0
    else if (ind%3==1)
    {
        S_curr_next=S_uni(coord_1_frame_curr[ind],coord_1_frame_next[ind],
            0,box_y_val,h);
        S_next_curr=S_uni(coord_1_frame_next[ind],coord_1_frame_curr[ind],
            0,box_y_val,h);
    }//end mod 3==1
    else if (ind %3==2)
    {
        S_curr_next=S_uni(coord_1_frame_curr[ind],coord_1_frame_next[ind],
            0,box_z_val,h);
        S_next_curr=S_uni(coord_1_frame_next[ind],coord_1_frame_curr[ind],
        0,box_z_val,h);

    }//end mod 3 ==2
    else
    {
        std::cerr<<"invalid ind="<<ind<<std::endl;
        std::exit(15);
    }

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


void mc_computation::energy_update_box_one_component(const std::vector<double>&box_1_frame_curr,
       const std::vector<double>&box_1_frame_next,const std::vector<double>&coord_1_frame_curr,double& UCurr, double& UNext)
{
    //energy curr
    this->dp_model.compute(UCurr,force_1_frame_curr,virial_1_frame_curr,coord_1_frame_curr,atom_type,box_1_frame_curr);

    //energy next
    this->dp_model.compute(UNext,force_1_frame_next,virial_1_frame_next,coord_1_frame_curr,atom_type,box_1_frame_next);
}


double mc_computation::box_2_volume(const std::vector<double>& box_1_frame)
{
    if (box_1_frame.size() !=box_component_num)
    {
        std::cerr<<"Incorrect box length: box_1_frame.size()="<<box_1_frame.size()<<std::endl;
        std::exit(16);
    }//end length check
    double box_x=box_1_frame[box_x_component_position];
    double box_y=box_1_frame[box_y_component_position];
    double box_z=box_1_frame[box_z_component_position];

    double volume=box_x*box_y*box_z;
    if (volume<=0)
    {
        std::cerr<<"invalid volume: volume="<<volume<<std::endl;
        std::exit(12);
    }//end if
    return volume;
}


void mc_computation::energy_update_coord_one_component(const std::vector<double>& coord_1_frame_curr,
        const std::vector<double>&coord_1_frame_next,  double& UCurr,   double& UNext,const std::vector<double>&box_1_frame_curr)
{
//energy curr
    this->dp_model.compute(UCurr,force_1_frame_curr,virial_1_frame_curr,coord_1_frame_curr,atom_type,box_1_frame_curr);

    //energy next
    this->dp_model.compute(UNext,force_1_frame_next,virial_1_frame_next,coord_1_frame_next,atom_type,box_1_frame_curr);




}




void  mc_computation::execute_mc_one_sweep(std::vector<double>& coord_1_frame_curr,
        std::vector<double>&box_1_frame_curr,
        double& UCurr,
       std::vector<double>&coord_1_frame_next,
       std::vector<double>&box_1_frame_next)
{
    double UNext = 0;
    double atom_x_max_val,atom_y_max_val,atom_z_max_val;
    // update coord
    for (int i=0;i<total_atom_xyz_components_num;i++)
    {
        int ind_coord=unif_in_0_NxNyNz_3_3_m1(e2);
        // std::cout<<"ind_coord="<<ind_coord<<std::endl;
        this->proposal_uni_coord(coord_1_frame_curr,coord_1_frame_next,ind_coord,box_1_frame_curr);
        this->energy_update_coord_one_component(coord_1_frame_curr,coord_1_frame_next,UCurr,UNext,box_1_frame_curr);

        double r_coord=this->acceptanceRatio_uni_for_coord(coord_1_frame_curr,coord_1_frame_next,
            ind_coord,UCurr,UNext,box_1_frame_curr[box_x_component_position],
            box_1_frame_curr[box_y_component_position],
            box_1_frame_curr[box_z_component_position]);
        double u = distUnif01(e2);
        if (u <= r_coord)
        {
            UCurr = UNext;
            std::memcpy(coord_1_frame_curr.data(),
                coord_1_frame_next.data(),total_atom_xyz_components_num*sizeof(double));

        }//end of accept-reject

    }//end updating coord
    //updating box
    for (int i=0;i<3;i++)
    {
        int which_ind_box=unif_in_0_2(e2);

        int ind_box=non_0_inds_in_box[which_ind_box];
        this->proposal_uni_box(box_1_frame_curr,box_1_frame_next,
            ind_box,
            coord_1_frame_curr,atom_x_max_val,
            atom_y_max_val,atom_z_max_val);

        this->energy_update_box_one_component(box_1_frame_curr,
            box_1_frame_next,coord_1_frame_curr,UCurr,UNext);

        double r_box=this->acceptanceRatio_uni_for_box(
            box_1_frame_curr,box_1_frame_next,ind_box,
            UCurr,UNext,
            atom_x_max_val,atom_y_max_val,atom_z_max_val);
        double u = distUnif01(e2);
        if (u<=r_box)
        {
            UCurr = UNext;
            std::memcpy(box_1_frame_curr.data(),
                box_1_frame_next.data(),box_component_num*sizeof
                (double));
        }//end of accept-reject


    }//end updating box
}