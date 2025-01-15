#importation
using JuMP
using Gurobi
include("lecture_distances.jl")
list_file = ["instance_6_1.txt"]
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
  model=Model(
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

  # On ne doit passer qu'une fois par aérodrome
  @constraint(model,[i in 1:n], sum(X[i,j] for j in (1:n)) <= 1 )

  # Distance maximale sur un arc est R 
  @constraint(model,[i in 1:n, j in 1:n], X[i,j]*D[i,j] <= R )

  # Nombre minimal d'aérodromes à visiter
  @constraint(model,sum(sum(X[i,j] for j in 1:n) for i in 1:n) >= Amin-1)

  # Contraintes sur les régions à visiter
  @constraint(model,[k in 1:Nr],sum(sum(X[i,j]+X[j,i] for i in 1:n) for j in regions[k]) >= 1)

  # Contrainte polynomiales
  list_i = deleteat!(collect(1:n), d)
  list_j = deleteat!(collect(1:n), d)
  @constraint(model,[i in list_i, j in list_j], u[j]>=u[i] + 1 - n*(1-X[i,j]))

  # Résolution
  JuMP.optimize!(model)
      
  # Affichage des résultats
  obj_value = JuMP.objective_value(model)
  println("Objective value: ", obj_value)
  println("Départ = ", d)
  println("Fin = ", f)

  for i in 1:n 
    for j in 1:n
      if JuMP.value(X[i,j]) == 1.0
        println("X[$i, $j] = 1 ")
      end
    end
  end
end