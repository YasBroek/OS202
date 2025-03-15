#include <string>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <thread>
#include <chrono>
#include <numeric>
#include <vector>
#include <omp.h>

#include "model.hpp"
#include "display.hpp"

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
        params.wind[0] = std::stoul(subkey);
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

int main( int nargs, char* args[] )
{
    auto params = parse_arguments(nargs-1, &args[1]);
    display_params(params);
    if (!check_params(params)) return EXIT_FAILURE;

    auto displayer = Displayer::init_instance( params.discretization, params.discretization );
    auto simu = Model( params.length, params.discretization, params.wind,
                       params.start);
    SDL_Event event;
    
    // Timing variables
    std::vector<double> update_times;
    std::vector<double> display_times;
    
    while (true)
    {
        // Measure model update time
        auto update_start = std::chrono::high_resolution_clock::now();
        bool updated = simu.update();
        auto update_end = std::chrono::high_resolution_clock::now();
        
        if (!updated) break;
        
        double update_duration = std::chrono::duration<double, std::milli>(update_end - update_start).count();
        update_times.push_back(update_duration);
        
        if ((simu.time_step() & 31) == 0) 
            std::cout << "Time step " << simu.time_step() << "\n===============" << std::endl;
        
        // Measure display update time
        auto display_start = std::chrono::high_resolution_clock::now();
        displayer->update( simu.vegetal_map(), simu.fire_map() );
        auto display_end = std::chrono::high_resolution_clock::now();
        
        double display_duration = std::chrono::duration<double, std::milli>(display_end - display_start).count();
        display_times.push_back(display_duration);
        
        // Print timing information periodically
        if ((simu.time_step() & 31) == 0) {
            double avg_update_time = std::accumulate(update_times.end() - std::min(static_cast<size_t>(32), update_times.size()), 
                                                    update_times.end(), 0.0) / 
                                    std::min(static_cast<size_t>(32), update_times.size());
            
            double avg_display_time = std::accumulate(display_times.end() - std::min(static_cast<size_t>(32), display_times.size()), 
                                                     display_times.end(), 0.0) / 
                                     std::min(static_cast<size_t>(32), display_times.size());
            
            std::cout << "Average update time: " << avg_update_time << " ms" << std::endl;
            std::cout << "Average display time: " << avg_display_time << " ms" << std::endl;
            std::cout << "Total step time: " << avg_update_time + avg_display_time << " ms" << std::endl;
        }
        
        if (SDL_PollEvent(&event) && event.type == SDL_QUIT)
            break;
        //std::this_thread::sleep_for(0.1s);
    }
    
    // Calculate and display final statistics
    if (!update_times.empty() && !display_times.empty()) {
        double total_update_time = std::accumulate(update_times.begin(), update_times.end(), 0.0);
        double total_display_time = std::accumulate(display_times.begin(), display_times.end(), 0.0);
        
        double avg_update_time = total_update_time / update_times.size();
        double avg_display_time = total_display_time / display_times.size();
        
        std::cout << "\n===== FINAL TIMING STATISTICS =====" << std::endl;
        std::cout << "Total time steps: " << update_times.size() << std::endl;
        std::cout << "Average model update time: " << avg_update_time << " ms" << std::endl;
        std::cout << "Average display update time: " << avg_display_time << " ms" << std::endl;
        std::cout << "Average total step time: " << avg_update_time + avg_display_time << " ms" << std::endl;
        std::cout << "Total simulation time: " << total_update_time + total_display_time << " ms" << std::endl;
    }
    
    return EXIT_SUCCESS;
}