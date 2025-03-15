#include <chrono>
#include <mpi.h>
#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <thread>
#include <chrono>

#include "model.hpp"
#include "display.hpp"

using namespace std::chrono;
using namespace std::string_literals;
using namespace std::chrono_literals;

struct ParamsType
{
    double length{1.};
    unsigned discretization{20u};
    std::array<double,2> wind{0.,0.};
    Model::LexicoIndices start{10u,10u};
};

void analyze_arg( int nargs, char* args[], ParamsType& params )
{
    if (nargs ==0) return;
    std::string key(args[0]);
    if (key == "-l"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une valeur pour la longueur du terrain !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.length = std::stoul(args[1]);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    auto pos = key.find("--longueur=");
    if (pos < key.size())
    {
        auto subkey = std::string(key,pos+11);
        params.length = std::stoul(subkey);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    if (key == "-n"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une valeur pour le nombre de cases par direction pour la discrétisation du terrain !" << std::endl;
            exit(EXIT_FAILURE);
        }
        params.discretization = std::stoul(args[1]);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--number_of_cases=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos+18);
        params.discretization = std::stoul(subkey);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    if (key == "-w"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une paire de valeurs pour la direction du vent !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values =std::string(args[1]);
        params.wind[0] = std::stod(values);
        auto pos = values.find(",");
        if (pos == values.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(values, pos+1);
        params.wind[1] = std::stod(second_value);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--wind=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos+7);
        params.wind[0] = std::stod(subkey);
        auto pos = subkey.find(",");
        if (pos == subkey.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(subkey, pos+1);
        params.wind[1] = std::stod(second_value);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }

    if (key == "-s"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Manque une paire de valeurs pour la position du foyer initial !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values =std::string(args[1]);
        params.start.column = std::stod(values);
        auto pos = values.find(",");
        if (pos == values.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la position du foyer initial" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(values, pos+1);
        params.start.row = std::stod(second_value);
        analyze_arg(nargs-2, &args[2], params);
        return;
    }
    pos = key.find("--start=");
    if (pos < key.size())
    {
        auto subkey = std::string(key, pos+8);
        params.start.column = std::stoul(subkey);
        auto pos = subkey.find(",");
        if (pos == subkey.size())
        {
            std::cerr << "Doit fournir deux valeurs séparées par une virgule pour définir la vitesse" << std::endl;
            exit(EXIT_FAILURE);
        }
        auto second_value = std::string(subkey, pos+1);
        params.start.row = std::stod(second_value);
        analyze_arg(nargs-1, &args[1], params);
        return;
    }
}

ParamsType parse_arguments( int nargs, char* args[] )
{
    if (nargs == 0) return {};
    if ( (std::string(args[0]) == "--help"s) || (std::string(args[0]) == "-h") )
    {
        std::cout << 
R"RAW(Usage : simulation [option(s)]
  Lance la simulation d'incendie en prenant en compte les [option(s)].
  Les options sont :
    -l, --longueur=LONGUEUR     Définit la taille LONGUEUR (réel en km) du carré représentant la carte de la végétation.
    -n, --number_of_cases=N     Nombre n de cases par direction pour la discrétisation
    -w, --wind=VX,VY            Définit le vecteur vitesse du vent (pas de vent par défaut).
    -s, --start=COL,ROW         Définit les indices I,J de la case où commence l'incendie (milieu de la carte par défaut)
)RAW";
        exit(EXIT_SUCCESS);
    }
    ParamsType params;
    analyze_arg(nargs, args, params);
    return params;
}

bool check_params(ParamsType& params)
{
    bool flag = true;
    if (params.length <= 0)
    {
        std::cerr << "[ERREUR FATALE] La longueur du terrain doit être positive et non nulle !" << std::endl;
        flag = false;
    }

    if (params.discretization <= 0)
    {
        std::cerr << "[ERREUR FATALE] Le nombre de cellules par direction doit être positive et non nulle !" << std::endl;
        flag = false;
    }

    if ( (params.start.row >= params.discretization) || (params.start.column >= params.discretization) )
    {
        std::cerr << "[ERREUR FATALE] Mauvais indices pour la position initiale du foyer" << std::endl;
        flag = false;
    }
    
    return flag;
}

void display_params(ParamsType const& params)
{
    std::cout << "Parametres définis pour la simulation : \n"
              << "\tTaille du terrain : " << params.length << std::endl 
              << "\tNombre de cellules par direction : " << params.discretization << std::endl 
              << "\tVecteur vitesse : [" << params.wind[0] << ", " << params.wind[1] << "]" << std::endl
              << "\tPosition initiale du foyer (col, ligne) : " << params.start.column << ", " << params.start.row << std::endl;
}

// Function for computation process
void computation_process(ParamsType& params) {
    // Initialize model
    auto simu = Model(params.length, params.discretization, params.wind, params.start);
    bool is_running = true;
    
    // Start time measurement
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Message buffer to send and receive commands/status
    int command;
    MPI_Status status;
    
    while (is_running) {
        // Wait for display process to request next update
        MPI_Recv(&command, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
        
        if (command == 0) { // Command to exit
            std::cout << "Computation process: Received exit command\n";
            break;
        }
        
        // Update simulation state
        is_running = simu.update();
        
        // Send simulation status
        int status_running = is_running ? 1 : 0;
        MPI_Send(&status_running, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        
        if (!is_running) {
            std::cout << "Computation process: Simulation naturally completed\n";
            // Espera um comando adicional do processo de exibição para garantir sincronização
            MPI_Recv(&command, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
            if (command == 0) {
                std::cout << "Computation process: Final confirmation received\n";
            }
            break;
        }
        
        // Send time step
        size_t time_step = simu.time_step();
        MPI_Send(&time_step, 1, MPI_UNSIGNED_LONG, 1, 0, MPI_COMM_WORLD);
        
        // Send vegetation map
        auto veg_map = simu.vegetal_map();
        MPI_Send(veg_map.data(), veg_map.size(), MPI_UINT8_T, 1, 0, MPI_COMM_WORLD);
        
        // Send fire map
        auto fire_map = simu.fire_map();
        MPI_Send(fire_map.data(), fire_map.size(), MPI_UINT8_T, 1, 0, MPI_COMM_WORLD);
        
        // Print progress info
        if ((time_step & 31) == 0) 
            std::cout << "Time step " << time_step << "\n===============" << std::endl;
    }
    
    // Record end time and print statistics
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end_time - start_time;
    std::cout << "Simulation time: " << duration.count() << " ms\n";
    std::cout << "Avg time per iteration: " 
              << (duration.count() / simu.time_step()) << " ms\n";
    
    std::cout << "Computation process completed successfully\n";
}

// Function for display process
void display_process(ParamsType& params) {
    // Initialize the displayer
    auto displayer = Displayer::init_instance(params.discretization, params.discretization);
    
    // Variables to receive simulation state
    bool is_running = true;
    size_t time_step;
    std::vector<std::uint8_t> vegetation_map(params.discretization * params.discretization);
    std::vector<std::uint8_t> fire_map(params.discretization * params.discretization);
    
    // For SDL events
    SDL_Event event;
    
    // Message buffer for commands/status
    int command;
    int status_running;
    MPI_Status mpi_status;
    
    while (true) {
        // Process SDL events
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                std::cout << "Display process: User closed window\n";
                // Signal computation process to exit
                command = 0;
                MPI_Send(&command, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                
                // Clean up resources
                displayer.reset();
                SDL_Quit();
                std::cout << "Display process completed successfully\n";
                return;
            }
        }
        
        // Request next simulation update
        command = 1;
        MPI_Send(&command, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        
        // Receive simulation status
        MPI_Recv(&status_running, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpi_status);
        is_running = (status_running == 1);
        
        if (!is_running) {
            // Simulation has ended
            std::cout << "Display process: Simulation has ended\n";
            
            // Pause for final visualization
            std::this_thread::sleep_for(2s);
            
            // Clean up resources
            displayer.reset();
            SDL_Quit();
            
            // Send final confirmation to computation process
            command = 0;
            MPI_Send(&command, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            
            std::cout << "Display process completed successfully\n";
            return;  // Sai da função display_process
        }
        
        // Receive time step
        MPI_Recv(&time_step, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD, &mpi_status);
        
        // Receive vegetation map
        MPI_Recv(vegetation_map.data(), vegetation_map.size(), MPI_UINT8_T, 
                0, 0, MPI_COMM_WORLD, &mpi_status);
        
        // Receive fire map
        MPI_Recv(fire_map.data(), fire_map.size(), MPI_UINT8_T, 
                0, 0, MPI_COMM_WORLD, &mpi_status);
        
        // Update display
        displayer->update(vegetation_map, fire_map);
    }
}

int main(int nargs, char* args[])
{
    // Initialize MPI
    int rank, size;
    MPI_Init(&nargs, &args);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // Check that we have exactly 2 processes
    if (size != 2) {
        if (rank == 0) {
            std::cerr << "This program requires exactly 2 MPI processes." << std::endl;
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    
    // Parse command line arguments
    auto params = parse_arguments(nargs-1, &args[1]);
    
    // Each process displays the parameters
    if (rank == 0) {
        std::cout << "Process " << rank << " (Computation): ";
        display_params(params);
    } else {
        std::cout << "Process " << rank << " (Display): ";
    display_params(params);
    }
    
    // Validate parameters
    if (!check_params(params)) {
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    
    try {
        // Branch based on process rank
        if (rank == 0) {
            computation_process(params);
        } else {
            display_process(params);
        }
    } catch (const std::exception& e) {
        std::cerr << "Process " << rank << " encountered an error: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Process " << rank << " encountered an unknown error" << std::endl;
    }
    
    // Finalize MPI
    std::cout << "Process " << rank << " finalizing MPI\n";
    MPI_Finalize();
    std::cout << "Process " << rank << " completed\n";
    return EXIT_SUCCESS;
}