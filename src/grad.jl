# Functions for preforming gradient descent on the likelihood of an embedding for a given graph

""" Likelihood of an edge given the embedding coordinate """
function edge_probability(f1, i::Int, j::Int, coordinates)
    yi, yj = coordinates[[i,j], :]
    return f1((norm(yi - yj))^2)
end

""" A likelihood function """
f1(x) = 1/(1+exp(x))

""" Derivative of the likelihood function """
f1prime(x) = -exp(x)/((exp(x)+1)^2)

""" Gradient of the log-likelihood """
function diff_log_sig(x,y)
    z = (norm(x - y))^2
    return (f1prime(z)/f1(z))*2*(x-y)
end

""" Gradient of the log-likelihood with negative edge samples """
function diff_f(x,y,neg_sample)
    x1 = diff_log_sig(x,y)
    x2 = diff_log_sig(y,x)
    for sample in neg_sample
        x2 = x2 + diff_log_sig(y,sample)
    end
    return [x1,x2]
end

"""
    diff_likelihood(coordinates,g,nodes, γ, weight2,num_negative_edges,negsource,negdestin) -> coordinates_revised
Gradient descent step on the embedding `coordinates` given the graph `g`.
"""
function diff_likelihood(coordinates,g,nodes, γ, weight2,num_negative_edges,negsource,negdestin)
    coordinates_revised = [[0.0,0.0] for x in nodes]
    for ed in edges(g)
        ei, ej = src(ed), dst(ed)

        indi = findall(x->x==ei, nodes)[1]
        indj = findall(x->x==ej, nodes)[1]

        yi = coordinates[indi, :]
        yj = coordinates[indj, :]

        x1 = diff_log_sig(yi,yj)
        x2 = diff_log_sig(yj,yi)

        neg_ind = findall(x->x==ej,negdestin)

        num_of_neg_ej = length(neg_ind)
        neg_ind_weight = weight2[neg_ind]
        neg_ind_weight_norm = neg_ind_weight./(sum(neg_ind_weight)+0.000001)

        for i in 1:num_negative_edges
            edge_ind = sample(neg_ind, Weights(neg_ind_weight_norm))

            indii = findall(x->x==negsource[edge_ind][1], nodes)[1]
            indjj = findall(x->x==negdestin[edge_ind][1], nodes)[1]

            yii = coordinates[indii, :]
            yjj = coordinates[indjj, :]

            x2 = x2 .+ diff_log_sig(yj,yii)
        end

        coordinates_revised[indi] = coordinates_revised[indi] .+ x1
        coordinates_revised[indj] = coordinates_revised[indj] .+ x2

    end

    coordinates_revised = [[x[1]/ne(g),x[2]/ne(g)] for x in coordinates_revised]

    return coordinates_revised
end
