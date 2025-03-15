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

// Nouvelle structure pour les informations de décomposition du domaine
struct DomainInfo {
    int rank;               // Rang du processus actuel
    int total_processes;    // Nombre total de processus de calcul
    unsigned global_rows;   // Nombre total de lignes dans la grille globale
    unsigned global_cols;   // Nombre total de colonnes dans la grille globale
    unsigned local_start_row; // Index de départ de la ligne locale (dans le domaine global)
    unsigned local_rows;    // Nombre de lignes dans le sous-domaine local
    unsigned local_cols;    // Nombre de colonnes dans le sous-domaine local (généralement égal à global_cols)
    
    // Calcule la ligne globale à partir de la ligne locale
    unsigned globalRow(unsigned local_row) const {
        return local_start_row + local_row;
    }
    
    // Vérifie si un index global est dans le domaine local
    bool containsGlobalRow(unsigned global_row) const {
        return global_row >= local_start_row && global_row < (local_start_row + local_rows);
    }
    
    // Convertit un index global en index local
    unsigned toLocalRow(unsigned global_row) const {
        return global_row - local_start_row;
    }
    
    // Initialise les informations du domaine pour un processus spécifique
    static DomainInfo create(int rank, int total_comp_processes, unsigned discretization) {
        DomainInfo info;
        info.rank = rank;
        info.total_processes = total_comp_processes;
        info.global_rows = discretization;
        info.global_cols = discretization;
        // Calculer le nombre de lignes par processus (division égale)
        unsigned base_rows = discretization / total_comp_processes;
        unsigned remainder = discretization % total_comp_processes;
        // Calculer l'index de départ et le nombre de lignes pour ce processus
        info.local_start_row = 0;
        for (int i = 0; i < rank; i++) {
            info.local_start_row += base_rows + (i < remainder ? 1 : 0);
        }
        // Ce processus reçoit une ligne supplémentaire si son rang est inférieur au reste
        info.local_rows = base_rows + (rank < remainder ? 1 : 0);
        info.local_cols = discretization;
        return info;
    }
};

void analyze_arg( int nargs, char* args[], ParamsType& params )
{
    if (nargs ==0) return;
    std::string key(args[0]);
    if (key == "-l"s)
    {
        if (nargs < 2)
        {
            std::cerr << "Valeur manquante pour la longueur du terrain !" << std::endl;
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
            std::cerr << "Valeur manquante pour le nombre de cases par direction pour la discrétisation du terrain !" << std::endl;
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
            std::cerr << "Paire de valeurs manquante pour la direction du vent !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values =std::string(args[1]);
        params.wind[0] = std::stod(values);
        auto pos = values.find(",");
        if (pos == values.size())
        {
            std::cerr << "Deux valeurs séparées par une virgule doivent être fournies pour définir la vitesse" << std::endl;
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
            std::cerr << "Deux valeurs séparées par une virgule doivent être fournies pour définir la vitesse" << std::endl;
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
            std::cerr << "Paire de valeurs manquante pour la position du foyer initial !" << std::endl;
            exit(EXIT_FAILURE);
        }
        std::string values =std::string(args[1]);
        params.start.column = std::stod(values);
        auto pos = values.find(",");
        if (pos == values.size())
        {
            std::cerr << "Deux valeurs séparées par une virgule doivent être fournies pour définir la position du foyer initial" << std::endl;
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
            std::cerr << "Deux valeurs séparées par une virgule doivent être fournies pour définir la vitesse" << std::endl;
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
R"RAW(Utilisation : simulation [option(s)]
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
        std::cerr << "[ERREUR FATALE] Le nombre de cellules par direction doit être positif et non nul !" << std::endl;
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
    std::cout << "Paramètres définis pour la simulation : \n"
              << "\tTaille du terrain : " << params.length << std::endl 
              << "\tNombre de cellules par direction : " << params.discretization << std::endl 
              << "\tVecteur vitesse : [" << params.wind[0] << ", " << params.wind[1] << "]" << std::endl
              << "\tPosition initiale du foyer (col, ligne) : " << params.start.column << ", " << params.start.row << std::endl;
}

// Version modifiée pour le calcul avec décomposition de domaine
void computation_process(ParamsType& params, int rank) {
    // Définir des constantes pour identifier les processus
    const int DISPLAY_RANK = 3;
    const int NUM_COMP_PROCESSES = 3;
    // Créer des informations sur le sous-domaine pour ce processus
    DomainInfo domain = DomainInfo::create(rank, NUM_COMP_PROCESSES, params.discretization);
    // Imprimer les informations du domaine
    std::cout << "Processus " << rank << " gérant les lignes " 
              << domain.local_start_row << " jusqu'à " 
              << (domain.local_start_row + domain.local_rows - 1) 
              << " (total: " << domain.local_rows << " lignes)" << std::endl;
    
    // Initialiser le sous-modèle local
    // Note : Ici, vous auriez besoin de modifier la classe Model pour supporter la décomposition de domaine
    // Pour l'instant, supposons que nous pouvons initialiser un modèle avec un sous-domaine spécifique
    auto simu = Model(params.length, params.discretization, params.wind, params.start);
    // Buffers pour l'échange de données entre processus
    std::vector<uint8_t> top_ghost_veg, bottom_ghost_veg;
    std::vector<uint8_t> top_ghost_fire, bottom_ghost_fire;
    // Taille du buffer d'une ligne
    const size_t line_size = domain.local_cols;
    // Allouer des buffers pour l'échange de données
    if (rank > 0) {  // Si ce n'est pas le premier processus, il a besoin d'une ligne fantôme supérieure
        top_ghost_veg.resize(line_size);
        top_ghost_fire.resize(line_size);
    }
    if (rank < NUM_COMP_PROCESSES - 1) {  // Si ce n'est pas le dernier processus, il a besoin d'une ligne fantôme inférieure
        bottom_ghost_veg.resize(line_size);
        bottom_ghost_fire.resize(line_size);
    }
    bool is_running = true;
    // Mesurer le temps d'exécution
    auto start_time = std::chrono::high_resolution_clock::now();
    // Buffers pour la communication
    int command;
    MPI_Status status;
    // Vecteurs pour stocker les cartes locales
    std::vector<uint8_t> local_veg_map(domain.local_rows * domain.local_cols);
    std::vector<uint8_t> local_fire_map(domain.local_rows * domain.local_cols);
    // Buffers pour collecter toutes les données pour le processus d'affichage
    std::vector<uint8_t> global_veg_map;
    std::vector<uint8_t> global_fire_map;
    // Si c'est le processus 0, allouez un buffer pour collecter toutes les données
    if (rank == 0) {
        global_veg_map.resize(params.discretization * params.discretization);
        global_fire_map.resize(params.discretization * params.discretization);
    }
    size_t iteration = 0;
    while (is_running) {
        // Tous les processus de calcul attendent le signal du processus d'affichage
        MPI_Recv(&command, 1, MPI_INT, DISPLAY_RANK, 0, MPI_COMM_WORLD, &status);
        if (command == 0) {
            std::cout << "Processus " << rank << ": Commande de sortie reçue\n";
            break;
        }
        // Incrémenter le compteur d'itérations
        iteration++;
        // Échanger les lignes fantômes entre les processus voisins
        // Seulement pour illustration, ici nous devrons adapter le modèle pour extraire ces données
        if (NUM_COMP_PROCESSES > 1) {
            // Envoyer la ligne inférieure au prochain processus
            if (rank < NUM_COMP_PROCESSES - 1) {
                // Obtenir la ligne inférieure actuelle du domaine local
                // Exemple : simu.get_bottom_row(bottom_ghost_veg, bottom_ghost_fire);
                // Extraire la dernière ligne de la carte locale
                for (size_t col = 0; col < domain.local_cols; col++) {
                    bottom_ghost_veg[col] = local_veg_map[(domain.local_rows - 1) * domain.local_cols + col];
                    bottom_ghost_fire[col] = local_fire_map[(domain.local_rows - 1) * domain.local_cols + col];
                }
                // Envoyer au prochain processus
                MPI_Send(bottom_ghost_veg.data(), line_size, MPI_UINT8_T, rank + 1, 10, MPI_COMM_WORLD);
                MPI_Send(bottom_ghost_fire.data(), line_size, MPI_UINT8_T, rank + 1, 11, MPI_COMM_WORLD);
            }
            // Recevoir la ligne supérieure du processus précédent
            if (rank > 0) {
                MPI_Recv(top_ghost_veg.data(), line_size, MPI_UINT8_T, rank - 1, 10, MPI_COMM_WORLD, &status);
                MPI_Recv(top_ghost_fire.data(), line_size, MPI_UINT8_T, rank - 1, 11, MPI_COMM_WORLD, &status);
                // Mettre à jour la ligne fantôme supérieure
                // Exemple : simu.set_top_ghost_row(top_ghost_veg, top_ghost_fire);
            }
            // Envoyer la ligne supérieure au processus précédent
            if (rank > 0) {
                // Obtenir la première ligne du domaine local
                for (size_t col = 0; col < domain.local_cols; col++) {
                    top_ghost_veg[col] = local_veg_map[col];
                    top_ghost_fire[col] = local_fire_map[col];
                }
                // Envoyer au processus précédent
                MPI_Send(top_ghost_veg.data(), line_size, MPI_UINT8_T, rank - 1, 20, MPI_COMM_WORLD);
                MPI_Send(top_ghost_fire.data(), line_size, MPI_UINT8_T, rank - 1, 21, MPI_COMM_WORLD);
            }
            // Recevoir la ligne inférieure du prochain processus
            if (rank < NUM_COMP_PROCESSES - 1) {
                MPI_Recv(bottom_ghost_veg.data(), line_size, MPI_UINT8_T, rank + 1, 20, MPI_COMM_WORLD, &status);
                MPI_Recv(bottom_ghost_fire.data(), line_size, MPI_UINT8_T, rank + 1, 21, MPI_COMM_WORLD, &status);
                // Mettre à jour la ligne fantôme inférieure
                // Exemple : simu.set_bottom_ghost_row(bottom_ghost_veg, bottom_ghost_fire);
            }
        }
    }
    // Mettre à jour la simulation locale
    bool local_running = true;
    // Ici, vous devez adapter le modèle pour travailler avec des sous-domaines
    // Exemple : local_running = simu.update_subdomain(domain, top_ghost_fire, bottom_ghost_fire);
    // Obtenir les cartes locales mises à jour
    // Exemple : 
    // simu.get_local_veg_map(local_veg_map);
    // simu.get_local_fire_map(local_fire_map);
    // Pour la simulation, supposons que nous calculons les nouveaux états
    // Normalement, cela serait fait par la classe Model modifiée
    // local_running = simu.update();
    // Rassembler l'état de tous les processus
    int global_status;
    int local_status = local_running ? 1 : 0;
    // Déterminer si la simulation doit continuer ou s'arrêter
    MPI_Allreduce(&local_status, &global_status, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    is_running = (global_status > 0);
    // Envoyer l'état au processus d'affichage
    if (rank == 0) {
        MPI_Send(&global_status, 1, MPI_INT, DISPLAY_RANK, 0, MPI_COMM_WORLD);
        if (!is_running) {
            std::cout << "Processus " << rank << ": Simulation terminée naturellement\n";
            // Attendre la confirmation du processus d'affichage
            MPI_Recv(&command, 1, MPI_INT, DISPLAY_RANK, 0, MPI_COMM_WORLD, &status);
            if (command == 0) {
                std::cout << "Processus " << rank << ": Confirmation finale reçue\n";
            }
        }
    }
    // Continuer uniquement si la simulation n'est pas terminée
    if (is_running) {
        // Envoyer le temps de simulation actuel (seulement le processus 0)
        if (rank == 0) {
            MPI_Send(&iteration, 1, MPI_UNSIGNED_LONG, DISPLAY_RANK, 0, MPI_COMM_WORLD);
        }
        // Collecter les données de tous les processus pour l'affichage
        // Chaque processus envoie ses données locales au processus 0
        if (rank == 0) {
            // Le processus 0 copie ses propres données locales
            for (size_t row = 0; row < domain.local_rows; row++) {
                for (size_t col = 0; col < domain.local_cols; col++) {
                    size_t local_idx = row * domain.local_cols + col;
                    size_t global_idx = (domain.local_start_row + row) * domain.global_cols + col;
                    global_veg_map[global_idx] = local_veg_map[local_idx];
                    global_fire_map[global_idx] = local_fire_map[local_idx];
                }
            }
            // Le processus 0 reçoit les données des autres processus
            for (int src_rank = 1; src_rank < NUM_COMP_PROCESSES; src_rank++) {
                // Obtenir les informations du domaine du processus source
                DomainInfo src_domain = DomainInfo::create(src_rank, NUM_COMP_PROCESSES, params.discretization);
                // Recevoir les données du processus source
                std::vector<uint8_t> recv_veg(src_domain.local_rows * src_domain.local_cols);
                std::vector<uint8_t> recv_fire(src_domain.local_rows * src_domain.local_cols);
                MPI_Recv(recv_veg.data(), recv_veg.size(), MPI_UINT8_T, src_rank, 30, MPI_COMM_WORLD, &status);
                MPI_Recv(recv_fire.data(), recv_fire.size(), MPI_UINT8_T, src_rank, 31, MPI_COMM_WORLD, &status);
                // Copier les données reçues dans les cartes globales
                for (size_t row = 0; row < src_domain.local_rows; row++) {
                    for (size_t col = 0; col < src_domain.local_cols; col++) {
                        size_t local_idx = row * src_domain.local_cols + col;
                        size_t global_idx = (src_domain.local_start_row + row) * domain.global_cols + col;
                        global_veg_map[global_idx] = recv_veg[local_idx];
                        global_fire_map[global_idx] = recv_fire[local_idx];
                    }
                }
            }
            // Envoyer les cartes globales au processus d'affichage
            MPI_Send(global_veg_map.data(), global_veg_map.size(), MPI_UINT8_T, DISPLAY_RANK, 0, MPI_COMM_WORLD);
            MPI_Send(global_fire_map.data(), global_fire_map.size(), MPI_UINT8_T, DISPLAY_RANK, 0, MPI_COMM_WORLD);
        } else {
            // Les autres processus envoient leurs données locales au processus 0
            MPI_Send(local_veg_map.data(), local_veg_map.size(), MPI_UINT8_T, 0, 30, MPI_COMM_WORLD);
            MPI_Send(local_fire_map.data(), local_fire_map.size(), MPI_UINT8_T, 0, 31, MPI_COMM_WORLD);
        }
        // Afficher les informations de progression
        if (rank == 0 && ((iteration & 31) == 0)) {
            std::cout << "Itération " << iteration << "\n===============" << std::endl;
        }
    }
}
// Enregistrer le temps final et imprimer les statistiques
auto end_time = std::chrono::high_resolution_clock::now();
std::chrono::duration<double, std::milli> duration = end_time - start_time;
std::cout << "Processus " << rank << ": Temps de simulation: " << duration.count() << " ms\n";
std::cout << "Processus " << rank << ": Temps moyen par itération: " 
          << (duration.count() / iteration) << " ms\n";
std::cout << "Processus " << rank << ": Calcul terminé avec succès\n";
}

// Fonction pour le processus d'affichage
void display_process(ParamsType& params) {
// Initialiser l'afficheur
auto displayer = Displayer::init_instance(params.discretization, params.discretization);
// Variables pour recevoir l'état de la simulation
bool is_running = true;
size_t time_step;
std::vector<std::uint8_t> vegetation_map(params.discretization * params.discretization);
std::vector<std::uint8_t> fire_map(params.discretization * params.discretization);
// Pour les événements SDL
SDL_Event event;
// Buffer pour les commandes/état
int command;
int status_running;
MPI_Status mpi_status;
while (true) {
    // Traiter les événements SDL
    while (SDL_PollEvent(&event)) {
        if (event.type == SDL_QUIT) {
            std::cout << "Processus d'affichage: L'utilisateur a fermé la fenêtre\n";
            // Signaler à tous les processus de calcul de sortir
            command = 0;
            for (int rank = 0; rank < 3; rank++) {
                MPI_Send(&command, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
            }
            // Nettoyer les ressources
            displayer.reset();
            SDL_Quit();
            std::cout << "Processus d'affichage terminé avec succès\n";
            return;
        }
    }
    // Demander la prochaine mise à jour de la simulation
    command = 1;
    for (int rank = 0; rank < 3; rank++) {
        MPI_Send(&command, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
    }
    // Recevoir l'état de la simulation (uniquement depuis le processus 0)
    MPI_Recv(&status_running, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpi_status);
    is_running = (status_running == 1);
    if (!is_running) {
        // Simulation terminée
        std::cout << "Processus d'affichage: Simulation terminée\n";
        // Pause pour la visualisation finale
        std::this_thread::sleep_for(2s);
        // Nettoyer les ressources
        displayer.reset();
        SDL_Quit();
        // Envoyer une confirmation finale à tous les processus de calcul
        command = 0;
        for (int rank = 0; rank < 3; rank++) {
            MPI_Send(&command, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
        }
        std::cout << "Processus d'affichage terminé avec succès\n";
        return;
    }
    // Recevoir le temps de simulation
    MPI_Recv(&time_step, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD, &mpi_status);
    // Recevoir les cartes de végétation et de feu du processus 0 (qui a déjà collecté de tous)
    MPI_Recv(vegetation_map.data(), vegetation_map.size(), MPI_UINT8_T, 
            0, 0, MPI_COMM_WORLD, &mpi_status);
    MPI_Recv(fire_map.data(), fire_map.size(), MPI_UINT8_T, 
            0, 0, MPI_COMM_WORLD, &mpi_status);
    // Mettre à jour l'affichage
    displayer->update(vegetation_map, fire_map);
}
}

int main(int nargs, char* args[])
{
// Initialiser MPI
int rank, size;
MPI_Init(&nargs, &args);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
// Vérifier que nous avons exactement 4 processus
if (size != 4) {
    if (rank == 0) {
        std::cerr << "Ce programme nécessite exactement 4 processus MPI." << std::endl;
    }
    MPI_Finalize();
    return EXIT_FAILURE;
}
// Analyser les arguments de la ligne de commande
auto params = parse_arguments(nargs-1, &args[1]);
// Chaque processus affiche les paramètres
if (rank < 3) {
    std::cout << "Processus " << rank << " (Calcul): ";
    display_params(params);
} else {
    std::cout << "Processus " << rank << " (Affichage): ";
    display_params(params);
}
// Valider les paramètres
if (!check_params(params)) {
    MPI_Finalize();
    return EXIT_FAILURE;
}
try {
    // Différencier entre les processus de calcul et d'affichage
    if (rank < 3) {
        computation_process(params, rank);
    } else {
        display_process(params);
    }
} catch (const std::exception& e) {
    std::cerr << "Processus " << rank << " a rencontré une erreur: " << e.what() << std::endl;
} catch (...) {
    std::cerr << "Processus " << rank << " a rencontré une erreur inconnue" << std::endl;
}
// Finaliser MPI
std::cout << "Processus " << rank << " finalise MPI\n";
MPI_Finalize();
std::cout << "Processus " << rank << " terminé\n";
return EXIT_SUCCESS;
}