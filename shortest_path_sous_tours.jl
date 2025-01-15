#importation
using JuMP
using Gurobi
include("lecture_distances.jl")

list_file = ["instance_20_1.txt"]
# ,"instance_20_1.txt", "instance_20_2.txt", "instance_20_3.txt", "instance_20_1.txt", 
# "instance_30_1.txt",
# "instance_40_1.txt", "instance_50_1.txt", "instance_70_1.txt", "instance_80_1.txt"]

for file_instance in list_file 
  #Get instance data
  n,d,f,Amin,Nr,R,regions,coords,D = readInstance(file_instance)

  V = zeros(n)
  V[d] = 1
  V[f] = -1

  # Déclaration du modèle
  model = Model(Gurobi.Optimizer)
  model = Model(
          optimizer_with_attributes(
            Gurobi.Optimizer, "Presolve" => 0
          ))
  # Déclaration des variables
  @variable(model, X[1:n,1:n] >=0, Bin)  #Arc Xij
  @variable(model, u[1:n]>=0, Int ) # u_i est nombre de sommets visités depuis le premier sommet
  # Déclaration de l'objectif
  @objective(model, Min, sum(sum(D[i,j]*X[i,j] for j in (1:n) for i in (1:n))))

  # Déclaration des contraintes :
  # Non-déplacement sur le même noeud
  @constraint(model,[t in 1:n], X[t,t] == 0)

  # Conservation du flot
  @constraint(model,[i in 1:n], sum(X[i,j]-X[j,i] for j in 1:n) == V[i])

  # On ne doit arriver qu'une fois à un aérodrome
  @constraint(model,[i in 1:n], sum(X[i,j] for j in (1:n)) <= 1 )
  # On ne doit partir qu'une fois à un aérodrome
  @constraint(model,[i in 1:n], sum(X[j,i] for j in (1:n)) <= 1 )
  # Distance maximale sur un arc est R 
  @constraint(model,[i in 1:n, j in 1:n], X[i,j]*D[i,j] <= R )

  # Nombre minimal d'aérodromes à visiter
  @constraint(model,sum(sum(X[i,j] for j in 1:n) for i in 1:n) >= Amin-1)

  # Contraintes sur les régions à visiter
  @constraint(model,[k in 1:Nr],sum(sum(X[i,j]+X[j,i] for i in 1:n) for j in regions[k]) >= 1)

  set_silent(model)

  # 1ère résolution
  JuMP.optimize!(model)

  global l = 0
  #initialisation du sous_problème 
  while  l == 0 
    X_tilde = JuMP.value.(X)
    # println("\nX_tilde = $X_tilde\n")

    # Déclaration du modèle de sous-tours
    model_st = Model(Gurobi.Optimizer)
    model_st =Model(
            optimizer_with_attributes(
            Gurobi.Optimizer, "Presolve" => 0
            ))
    set_silent(model_st)
    
    # Déclaration des variables
    @variable(model_st, a[1:n] >=0, Bin)  
    # Déclaration de la fonction objective
    @objective(model_st, Max, sum(sum(X_tilde[i,j]*a[i]*a[j]  for j in 1:n )  for i in 1:n ) - sum(a[i] for i in 1:n) + 1 )
    # Déclaration des contraintes
    @constraint(model_st,sum(a[i] for i in 1:n ) >= 2 )
    @constraint(model_st,sum(a[i] for i in 1:n) <= n-2 )

    println("\n\n========================= Résolution du problème sous-tours ==================================\n")
    # Optimisation
    # println(model_st)
    JuMP.optimize!(model_st)
    println(("solution sous-tour = $(JuMP.objective_value(model_st)) ") )
    if JuMP.objective_value(model_st) <= 0.5 # ie pas de sous-tours trouvés
      global l = 1
    else 
      global l = 0
      a_tilde = JuMP.value.(a)
      println("a_tilde (sous-tour trouvé) = ", a_tilde)
      # Ajout nouvelle contrainte au modèle de base
      for k in 1:n 
        if a_tilde[k] == 1 
          list_i = deleteat!(collect(1:n), k)
          @constraint(model,sum(sum(X[i,j]*a_tilde[i]*a_tilde[j] for j in 1:n ) for i in 1:n) <= sum(sum(X[i,j]*a_tilde[i] for j in 1:n) for i in list_i))
        end
      end
      println("\n\n========================= Résolution du nouveau problème ==================================\n")
      JuMP.optimize!(model)
      
    end
  end 


  # Affichage des résultats finaux
  obj_value = JuMP.objective_value(model)
  println("Objective value: ", obj_value)
  println("Départ = ", d)
  println("Fin = ", f)
  println("Node count = ", JuMP.node_count(model))
  println("Solve time = ", JuMP.solve_time(model))

  for i in 1:n 
    for j in 1:n
      if JuMP.value(X[i,j]) == 1.0
        println("X[$i, $j] = 1 ")
      end
    end
  end

  # Ecriture des résultats dans 1 fichier .txt
  open("report_exponential.txt", "a") do file
    write(file, "$file_instance \n") # instance
    write(file, "Départ = $d \n")
    write(file, "Objective value: $obj_value")
    write(file, "Fin = $f\n")
    write(file, "Amin = $Amin\n")
    write(file, "Node count = $node_count\n")
    write(file, "Solve time = $solve_time\n")
    # write(file, "Simplex iterations = ", string(JuMP.simplex_iterations(model)))
    # write(file, summary)
    write(file, "================ END ====================\n \n \n")
  end

end