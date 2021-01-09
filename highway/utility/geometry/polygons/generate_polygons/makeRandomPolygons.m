function P = makeRandomPolygons(Npolygons)
    Nvertices = round(randRange(3,9,[],[],1,Npolygons)) ;
    P = makeRandomPolygon(Nvertices) ;
end