import math
import random
import time
import matplotlib.pyplot as plt

# =============================================================================
# SECTION 1 : UTILITAIRES ARITHMÉTIQUES (MPI - RIGOUREUX)
# =============================================================================

def algorithme_euclide_etendu(a, b):
    """
    Calcule le PGCD et les coefficients de Bézout.
    ax + by = pgcd(a, b)
    """
    if a == 0:
        return b, 0, 1
    pgcd, x1, y1 = algorithme_euclide_etendu(b % a, a)
    x = y1 - (b // a) * x1
    y = x1
    return pgcd, x, y

def inverse_modulaire(a, m):
    """
    Calcule l'inverse de a modulo m via Euclide étendu.
    """
    pgcd, x, y = algorithme_euclide_etendu(a, m)
    if pgcd != 1:
        return None
    return x % m

def resoudre_congruence_lineaire(v, u, q):
    """
    Résout l'équation v * x = u (mod q).
    Utilisé pour extraire le secret après la collision.
    """
    pgcd, x0, y0 = algorithme_euclide_etendu(v, q)
    if u % pgcd != 0:
        return [] # Aucune solution possible
    
    solutions = []
    # Solution particulière mod (q/pgcd)
    base_sol = (x0 * (u // pgcd)) % (q // pgcd)
    # Il y a exactement 'pgcd' solutions dans Z/qZ
    for i in range(pgcd):
        solutions.append((base_sol + i * (q // pgcd)) % q)
    return solutions

def est_premier(n):
    """ Test de primalité déterministe pour les petites/moyennes tailles. """
    if n < 2: return False
    if n == 2 or n == 3: return True
    if n % 2 == 0: return False
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        if n % i == 0: return False
    return True

# =============================================================================
# SECTION 2 : MARCHE ALÉATOIRE DE TESKE (3 PARTITIONS)
# =============================================================================

def etape_teske(x, a, b, g, h, p, q):
    """
    Réalise une étape de la marche de Teske.
    Invariant maintenu : x = g^a * h^b (mod p).
    """
    # La zone est déterminée par le reste de la division par 3
    zone = x % 3
    
    if zone == 0:
        # Zone S0 : Élévation au carré (Doublage des exposants)
        nouveau_x = (x * x) % p
        nouveau_a = (a * 2) % q
        nouveau_b = (b * 2) % q
    elif zone == 1:
        # Zone S1 : Multiplication par le générateur g
        nouveau_x = (x * g) % p
        nouveau_a = (a + 1) % q
        nouveau_b = b % q
    else:
        # Zone S2 : Multiplication par la cible h
        nouveau_x = (x * h) % p
        nouveau_a = a % q
        nouveau_b = (b + 1) % q
        
    return nouveau_x, nouveau_a, nouveau_b

# =============================================================================
# SECTION 3 : MÉTHODE EXHAUSTIVE (BRUTE FORCE - O(q))
# =============================================================================

def resoudre_dlp_exhaustif(g, h, p, q):
    """
    Recherche exhaustive du secret x en testant chaque possibilité de 0 à q-1.
    Complexité temporelle : O(q).
    """
    nb_evaluations = 0
    # On parcourt tout l'espace des clés possibles
    for x_candidat in range(q):
        nb_evaluations += 1
        # Test de validité : g^x == h (mod p)
        if pow(g, x_candidat, p) == h:
            return x_candidat, nb_evaluations
    return None, nb_evaluations



# =============================================================================
# SECTION 3 : RÉSOLUTION DU DLP - MÉTHODE DE FLOYD
# =============================================================================

def resoudre_dlp_floyd(g, h, p, q):
    """
    Résout g^x = h (mod p) via l'algorithme de Floyd.
    Complexité théorique : ~3.09 * sqrt(q) évaluations de f.
    """
    # Initialisation : x = g^0 * h^0 = 1
    xt, at, bt = 1, 0, 0 # Tortue
    xl, al, bl = 1, 0, 0 # Lièvre
    
    nb_evaluations = 0
    
    while True:
        # La tortue avance de 1 pas
        xt, at, bt = etape_teske(xt, at, bt, g, h, p, q)
        # Le lièvre avance de 2 pas
        xl, al, bl = etape_teske(xl, al, bl, g, h, p, q)
        xl, al, bl = etape_teske(xl, al, bl, g, h, p, q)
        
        nb_evaluations += 3 # 1 pour T, 2 pour L
        
        if xt == xl:
            # Collision détectée !
            v = (bl - bt) % q
            u = (at - al) % q
            
            solutions = resoudre_congruence_lineaire(v, u, q)
            for x_cand in solutions:
                if pow(g, x_cand, p) == h:
                    return x_cand, nb_evaluations
            # Si échec de résolution (rare avec Teske), on continue la boucle
            # pour trouver une autre collision ou on pourrait réinitialiser x0

# =============================================================================
# SECTION 4 : RÉSOLUTION DU DLP - MÉTHODE DE BRENT
# =============================================================================

def resoudre_dlp_brent(g, h, p, q):
    """
    Résout g^x = h (mod p) via l'algorithme de Brent.
    Complexité théorique : ~1.98 * sqrt(q) évaluations de f.
    """
    x_fixe, a_fixe, b_fixe = 1, 0, 0
    x_actu, a_actu, b_actu = 1, 0, 0
    
    puissance = 1
    nb_evaluations = 0
    
    while True:
        # Mise à jour du point fixe (balise) à chaque puissance de 2
        x_fixe, a_fixe, b_fixe = x_actu, a_actu, b_actu
        
        # On explore la fenêtre de taille 2^k
        for i in range(puissance):
            x_actu, a_actu, b_actu = etape_teske(x_actu, a_actu, b_actu, g, h, p, q)
            nb_evaluations += 1
            
            if x_actu == x_fixe:
                # Collision détectée !
                v = (b_actu - b_fixe) % q
                u = (a_fixe - a_actu) % q
                
                solutions = resoudre_congruence_lineaire(v, u, q)
                for x_cand in solutions:
                    if pow(g, x_cand, p) == h:
                        return x_cand, nb_evaluations
        
        puissance *= 2

# =============================================================================
# SECTION 5 : VALIDATION STATISTIQUE (FLAJOLET-ODLYZKO)
# =============================================================================

def verifier_statistiques_flajolet(nb_bits, nb_essais=5000):
    """
    Vérifie que la marche de Teske produit des structures de Rho conformes
    aux espérances de Flajolet et Odlyzko : E(mu) = E(lambda) = sqrt(pi*N/8).
    """
    print(f"\n[STATISTIQUES] Analyse de la structure en Rho sur {nb_bits} bits...")
    
    p = 0
    while not est_premier(p) or p < 100:
        p = random.getrandbits(nb_bits)
    q = p - 1
    g = 2
    h = pow(g, random.randint(2, q-1), p)
    
    queues = []
    cycles = []
    
    for _ in range(nb_essais):
        historique = {}
        x, a, b = 1, 0, 0
        i = 0
        while x not in historique:
            historique[x] = i
            x, a, b = etape_teske(x, a, b, g, h, p, q)
            i += 1
        
        mu = historique[x] # Longueur de la queue
        lam = i - historique[x] # Longueur du cycle
        queues.append(mu)
        cycles.append(lam)
        
    avg_mu = sum(queues) / nb_essais
    avg_lam = sum(cycles) / nb_essais
    theorie = math.sqrt((math.pi * q) / 8)
    
    print(f"  > Moyenne Queue (Exp) : {avg_mu:.2f} | Théorie : {theorie:.2f}")
    print(f"  > Moyenne Cycle (Exp) : {avg_lam:.2f} | Théorie : {theorie:.2f}")
    ecart = (abs(avg_mu + avg_lam - 2*theorie) / theorie) * 100
    print(f"  > Conformité statistique : {100 - ecart:.2f}%")

# =============================================================================
# SECTION 6 : BENCHMARK EXHAUSTIF ET MESURE DU GAIN DE BRENT
# =============================================================================

def executer_comparaison_performance(essais = 15):
    
    tailles_bits = [14, 16, 18, 20, 22]
    
    stats_floyd = []
    stats_brent = []
    gains_pourcentage = []
    
    print("\n" + "="*60)
    print(f"{'BITS':<6} | {'ÉVALS FLOYD':<15} | {'ÉVALS BRENT':<15} | {'GAIN (%)':<10}")
    print("-" * 60)
    
    for b in tailles_bits:
        total_eval_f = 0
        total_eval_b = 0
        
        for _ in range(essais):
            # Génération d'un problème DLP
            p = 0
            while not est_premier(p):
                p = random.getrandbits(b)
            q = p - 1
            g = random.randint(2, p-2)
            secret = random.randint(2, q-1)
            h = pow(g, secret, p)
            
            # Résolution Floyd
            _, ev_f = resoudre_dlp_floyd(g, h, p, q)
            total_eval_f += ev_f
            
            # Résolution Brent
            _, ev_b = resoudre_dlp_brent(g, h, p, q)
            total_eval_b += ev_b
            
        moy_f = total_eval_f / essais
        moy_b = total_eval_b / essais
        
        # Calcul du gain chiffré : (Floyd - Brent) / Floyd
        gain = ((moy_f - moy_b) / moy_f) * 100
        
        stats_floyd.append(moy_f)
        stats_brent.append(moy_b)
        gains_pourcentage.append(gain)
        
        print(f"{b:<6} | {moy_f:<15.1f} | {moy_b:<15.1f} | {gain:>8.2f}%")

    print("-" * 60)
    gain_moyen_total = sum(gains_pourcentage) / len(gains_pourcentage)
    print(f"GAIN D'EFFICACITÉ MOYEN DE BRENT : {gain_moyen_total:.2f}%")
    print("="*60)

    # --- GÉNÉRATION DU GRAPHIQUE ---
    plt.figure(figsize=(10, 6))
    #plt.plot(tailles_bits, 
    plt.plot(tailles_bits, stats_floyd, 'ro-', label='Floyd (Tortue/Lièvre)')
    plt.plot(tailles_bits, stats_brent, 'bs-', label='Brent (Fenêtrage)')
    plt.yscale('log')
    plt.xlabel('Taille du groupe (Bits)')
    plt.ylabel('Nombre moyen d\'évaluations de la fonction de marche (log)')
    plt.title('Démonstration de l\'efficacité supérieure de l\'algorithme de Brent')
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.5)
    
    # Ajout d'une zone mettant en évidence le gain
    plt.text(tailles_bits[1], stats_floyd[-1], f"Gain moyen observé : {gain_moyen_total:.1f}%", 
             bbox=dict(facecolor='green', alpha=0.2), fontsize=12)
    
    plt.show()
    
def executer_benchmark_temporel(essais = 15):
    """ 
    Mesure le temps d'exécution réel (CPU Time) pour mettre en évidence 
    la différence de complexité O(q) vs O(sqrt(q)).
    """
    # Tailles binaires pour le test (limité pour l'exhaustif)
    tailles_bits = [14, 16, 18, 20, 22] 
    
    temps_exh = []
    temps_floyd = []
    temps_brent = []
    
    print("\n" + "="*85)
    print(f"{'BITS':<6} | {'TEMPS EXHAUSTIF (s)':<20} | {'TEMPS FLOYD (s)':<18} | {'TEMPS BRENT (s)'}")
    print("-" * 85)
    
    for b in tailles_bits:
        t_exh, t_floyd, t_brent = 0.0, 0.0, 0.0
        for _ in range(essais):
            # Génération d'un problème DLP aléatoire
            p = 0
            while not est_premier(p): p = random.getrandbits(b)
            q = p - 1
            g = random.randint(2, p-2)
            h = pow(g, random.randint(2, q-1), p)
            
            # 1. Mesure Exhaustif
            start = time.perf_counter()
            resoudre_dlp_exhaustif(g, h, p, q)
            t_exh += (time.perf_counter() - start)
            
            # 2. Mesure Floyd
            start = time.perf_counter()
            resoudre_dlp_floyd(g, h, p, q)
            t_floyd += (time.perf_counter() - start)
            
            # 3. Mesure Brent
            start = time.perf_counter()
            resoudre_dlp_brent(g, h, p, q)
            t_brent += (time.perf_counter() - start)
            
        moy_e, moy_f, moy_b = t_exh/essais, t_floyd/essais, t_brent/essais
        temps_exh.append(moy_e)
        temps_floyd.append(moy_f)
        temps_brent.append(moy_b)
        
        print(f"{b:<6} | {moy_e:<20.6f} | {moy_f:<18.6f} | {moy_b:<15.6f}")

    # --- GÉNÉRATION DU GRAPHIQUE TEMPOREL ---
    plt.figure(figsize=(10, 6))
    plt.plot(tailles_bits, temps_exh, 'ko--', label='Recherche Exhaustive O(q)')
    plt.plot(tailles_bits, temps_floyd, 'ro-', label='Pollard Rho (Floyd) O(√q)')
    plt.plot(tailles_bits, temps_brent, 'bs-', label='Pollard Rho (Brent) O(√q)')
    
    plt.yscale('log') # Indispensable pour voir les pentes
    plt.xlabel('Taille du secret (Bits)')
    plt.ylabel('Temps d\'exécution moyen (secondes) - Échelle Log')
    plt.title('Analyse de la Complexité Temporelle Réelle')
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.5)
    
    
    plt.show()
    
def comparaison_finale_exhaustive_vs_pollard():
    """
    Compare l'approche exhaustive et Pollard Rho sur une instance DH de 24 bits.
    Affiche le facteur d'accélération pour la conclusion du TIPE.
    """
    print("\n" + "="*70)
    print(" COMPARISON : FORCE BRUTE vs POLLARD RHO (24 BITS)")
    print("="*70)

    # 1. Génération de l'instance
    bits = 24
    p = 0
    while not est_premier(p):
        p = random.randint(2**(bits-1), 2**bits - 1)
    q = p - 1
    g = 2
    secret = random.randint(2, q - 1)
    h = pow(g, secret, p)

    print(f"Instance : 2^x = {h} [mod {p}]")

    # 2. Résolution par Force Brute (Exhaustive)
    print("\n[1/2] Lancement Force Brute...")
    start_exh = time.perf_counter()
    # On simule la boucle de recherche pour compter les itérations
    # (resoudre_dlp_exhaustif retourne le secret et le nb d'itérations)
    sol_exh, it_exh = resoudre_dlp_exhaustif(g, h, p, q)
    end_exh = time.perf_counter()
    t_exh = end_exh - start_exh

    # 3. Résolution par Pollard Rho (Brent)
    print("[2/2] Lancement Pollard Rho (Brent)...")
    start_pollard = time.perf_counter()
    sol_p, it_pollard = resoudre_dlp_brent(g, h, p, q)
    end_pollard = time.perf_counter()
    t_pollard = end_pollard - start_pollard

    # 4. Calcul du gain
    facteur_vitesse = t_exh / t_pollard if t_pollard > 0 else 0
    reduction_operations = it_exh / it_pollard

    print("\n" + "-"*30 + " BILAN " + "-"*30)
    print(f"{'MÉTHODE':<20} | {'ITÉRATIONS':<15} | {'TEMPS (s)':<10}")
    print(f"{'Force Brute':<20} | {it_exh:<15,} | {t_exh:.4f}s")
    print(f"{'Pollard Rho':<20} | {it_pollard:<15,} | {t_pollard:.4f}s")
    print("-" * 70)
    print(f" ACCÉLÉRATION : x{int(facteur_vitesse)} fois plus rapide")
    print(f" RÉDUCTION DES CALCULS : facteur {int(reduction_operations)}")
    print(f" RÉSULTAT : x = {int(sol_p)}")
    print("="*70)

# =============================================================================
# POINT D'ENTRÉE DU PROGRAMME
# =============================================================================

if __name__ == "__main__":
    print("PROJET TIPE : ÉTUDE ET OPTIMISATION DE LA DÉTECTION DE CYCLE")
    print("Auteur : Raphaël DEPUYDT")
    
    # 1. Vérification de la structure du Rho (Flajolet)
    verifier_statistiques_flajolet(12, nb_essais=300)
    
    # 2. Benchmark et démonstration du gain de 36%
    executer_comparaison_performance(300)

    executer_benchmark_temporel()
    comparaison_finale_exhaustive_vs_pollard()