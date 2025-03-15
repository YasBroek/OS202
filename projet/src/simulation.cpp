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

// Processo de computação
void computation_process(ParamsType& params) {
    // Inicializa o modelo
    auto simu = Model(params.length, params.discretization, params.wind, params.start);
    bool is_running = true;

    // Inicia a medição de tempo
    auto start_time = high_resolution_clock::now();
    
    // Track model update times
    double total_update_time = 0.0;
    size_t update_count = 0;
    
    // Get number of hardware threads available
    unsigned int num_threads = std::thread::hardware_concurrency();
    std::cout << "Running with " << num_threads << " hardware threads available" << std::endl;

    // Buffer para enviar/receber comandos/status
    int command;
    MPI_Status status;

    // Variáveis para armazenar estado futuro
    bool next_state_ready = false;
    size_t next_time_step = 0;
    std::vector<std::uint8_t> next_veg_map(params.discretization * params.discretization);
    std::vector<std::uint8_t> next_fire_map(params.discretization * params.discretization);
    
    // Variables to track time advancement per thread
    double total_simulation_time = 0.0;
    double thread_advancement_rate = 0.0;

    while (is_running) {
        // Calcula o próximo estado em segundo plano
        if (!next_state_ready) {
            auto update_start = high_resolution_clock::now();
            next_state_ready = simu.update();
            auto update_end = high_resolution_clock::now();
            
            // Calculate and accumulate the update time
            std::chrono::duration<double, std::milli> update_duration = update_end - update_start;
            total_update_time += update_duration.count();
            update_count++;
            
            next_time_step = simu.time_step();
            next_veg_map = simu.vegetal_map();
            next_fire_map = simu.fire_map();
        }

        // Espera pelo comando do processo de exibição
        MPI_Recv(&command, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);

        if (command == 0) { // Comando para sair
            std::cout << "Computation process: Received exit command\n";
            break;
        }

        // Envia o estado atual para o display
        int status_running = next_state_ready ? 1 : 0;
        MPI_Send(&status_running, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);

        if (!next_state_ready) {
            std::cout << "Computation process: Simulation naturally completed\n";
            MPI_Recv(&command, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
            if (command == 0) {
                std::cout << "Computation process: Final confirmation received\n";
            }
            break;
        }

        // Envia o estado atual para o display
        MPI_Send(&next_time_step, 1, MPI_UNSIGNED_LONG, 1, 0, MPI_COMM_WORLD);
        MPI_Send(next_veg_map.data(), next_veg_map.size(), MPI_UINT8_T, 1, 0, MPI_COMM_WORLD);
        MPI_Send(next_fire_map.data(), next_fire_map.size(), MPI_UINT8_T, 1, 0, MPI_COMM_WORLD);

        // Libera o estado futuro para cálculo
        next_state_ready = false;

        // Imprime informações de progresso
        if ((next_time_step & 31) == 0)
            std::cout << "Time step " << next_time_step << "\n===============" << std::endl;
    }

    // Registra o tempo final e imprime estatísticas
    auto end_time = high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end_time - start_time;
    total_simulation_time = duration.count();
    
    std::cout << "Simulation time: " << total_simulation_time << " ms\n";
    std::cout << "Avg time per iteration: " << (total_simulation_time / simu.time_step()) << " ms\n";
    
    // Display the average model update time
    if (update_count > 0) {
        double avg_update_time = total_update_time / update_count;
        std::cout << "Avg model update time: " << avg_update_time << " ms\n";
        
        // Calculate thread advancement metrics
        if (num_threads > 0) {
            thread_advancement_rate = (simu.time_step() / total_simulation_time) * 1000.0; // time steps per second
            double per_thread_advancement = thread_advancement_rate / num_threads;
            std::cout << "Time advancement rate: " << thread_advancement_rate << " steps/second\n";
            std::cout << "Time advancement per thread: " << per_thread_advancement << " steps/second/thread\n";
            
            // Calculate theoretical speedup efficiency
            double ideal_time = total_update_time / num_threads;
            double speedup = total_update_time / avg_update_time;
            double efficiency = speedup / num_threads * 100.0;
            std::cout << "Parallel efficiency: " << efficiency << "%\n";
            
            // Interpret the results
            std::cout << "\nPerformance Analysis:\n";
            std::cout << "====================\n";
            if (efficiency > 90.0) {
                std::cout << "Excellent scaling efficiency. The simulation is taking full advantage of all threads.\n";
            } else if (efficiency > 70.0) {
                std::cout << "Good scaling efficiency. The simulation is benefiting well from parallelization.\n";
            } else if (efficiency > 50.0) {
                std::cout << "Moderate scaling efficiency. Some parallel overhead or load imbalance may be present.\n";
            } else {
                std::cout << "Poor scaling efficiency. The simulation may be limited by sequential bottlenecks,\n";
                std::cout << "communication overhead, or memory bandwidth constraints.\n";
            }
            
            // Analyze the relationship between model update time and overall simulation time
            double update_portion = total_update_time / total_simulation_time * 100.0;
            std::cout << "\nModel updates account for " << update_portion << "% of total simulation time.\n";
            if (update_portion < 50.0) {
                std::cout << "Most time is spent on communication or display. Improving model parallelization\n";
                std::cout << "may have limited impact on overall performance.\n";
            } else {
                std::cout << "Model computation dominates simulation time. Optimizing model parallelization\n";
                std::cout << "should yield significant overall performance improvements.\n";
            }
        }
    } else {
        std::cout << "No model updates performed\n";
    }
    
    std::cout << "Computation process completed successfully\n";
}

// Processo de exibição
void display_process(ParamsType& params) {
    // Inicializa o displayer
    auto displayer = Displayer::init_instance(params.discretization, params.discretization);

    // Variáveis para receber o estado da simulação
    bool is_running = true;
    size_t time_step;
    std::vector<std::uint8_t> vegetation_map(params.discretization * params.discretization);
    std::vector<std::uint8_t> fire_map(params.discretization * params.discretization);

    // Para eventos SDL
    SDL_Event event;

    // Buffer para comandos/status
    int command;
    int status_running;
    MPI_Status mpi_status;

    while (true) {
        // Processa eventos SDL
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                std::cout << "Display process: User closed window\n";
                command = 0;
                MPI_Send(&command, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
                displayer.reset();
                SDL_Quit();
                std::cout << "Display process completed successfully\n";
                return;
            }
        }

        // Solicita a próxima atualização da simulação
        command = 1;
        MPI_Send(&command, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

        // Recebe o status da simulação
        MPI_Recv(&status_running, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpi_status);
        is_running = (status_running == 1);

        if (!is_running) {
            std::cout << "Display process: Simulation has ended\n";
            std::this_thread::sleep_for(2s);
            displayer.reset();
            SDL_Quit();
            command = 0;
            MPI_Send(&command, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            std::cout << "Display process completed successfully\n";
            return;
        }

        // Recebe o passo de tempo
        MPI_Recv(&time_step, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD, &mpi_status);

        // Recebe o mapa de vegetação
        MPI_Recv(vegetation_map.data(), vegetation_map.size(), MPI_UINT8_T, 0, 0, MPI_COMM_WORLD, &mpi_status);

        // Recebe o mapa de fogo
        MPI_Recv(fire_map.data(), fire_map.size(), MPI_UINT8_T, 0, 0, MPI_COMM_WORLD, &mpi_status);

        // Atualiza a exibição
        displayer->update(vegetation_map, fire_map);
    }
}

int main(int nargs, char* args[]) {
    // Inicializa MPI
    int rank, size;
    MPI_Init(&nargs, &args);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Verifica se há exatamente 2 processos
    if (size != 2) {
        if (rank == 0) {
            std::cerr << "This program requires exactly 2 MPI processes.\n";
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    // Parseia os argumentos da linha de comando
    auto params = parse_arguments(nargs - 1, &args[1]);

    // Cada processo exibe os parâmetros
    if (rank == 0) {
        std::cout << "Process " << rank << " (Computation): ";
        display_params(params);
    } else {
        std::cout << "Process " << rank << " (Display): ";
        display_params(params);
    }

    // Valida os parâmetros
    if (!check_params(params)) {
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    // Executa o processo correspondente
    if (rank == 0) {
        computation_process(params);
    } else {
        display_process(params);
    }

    // Finaliza MPI
    std::cout << "Process " << rank << " finalizing MPI\n";
    MPI_Finalize();
    std::cout << "Process " << rank << " completed\n";

    return EXIT_SUCCESS;
}