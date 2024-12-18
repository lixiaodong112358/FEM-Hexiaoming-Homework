function result=generate_initial_vector(function_initial,Pb)

result=feval(function_initial,Pb(1,:),Pb(2,:))';


end