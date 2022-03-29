function [f,v] = field2mesh(Q,Params)

%set up meshgrid of points
[X,Y,Z] = meshgrid(Params.Extents(1,1):Params.dx:Params.Extents(1,2),...
    Params.Extents(2,1):Params.dx:Params.Extents(2,2),...
    Params.Extents(3,1):Params.dx:Params.Extents(3,2));

if Params.CapFaces ~= 0
    Q([1 end],:,:) = Params.CapFaces;
    Q(:,[1 end],:) = Params.CapFaces;
    Q(:,:,[1 end]) = Params.CapFaces;
end

%extract isosurface through field
p = isosurface(X,Y,Z,Q,Params.IsoVal);

f = p.faces; v = p.vertices;

valid_mesh = false;
fail = false;
if ~isempty(f) && ~isempty(v) && Params.Remesh

    %check the number of connected components
    if Params.rmIslands
        [NR,RS] = segment_connected_components(f, 'explicit');
    else
        NR = 1;
    end
    %if there are multiple regions, quit
    if NR==1
        if Params.CloseHoles
            %close any holes in the mesh
            try
                [f,v] = triSurfCloseHoles(f,v);
            catch
                fprintf('\nClose Holes Failed\n')
                fail = true;
            end
        end

        if Params.Remesh && ~fail
            %remesh the isosurface
            opt.pointSpacing=Params.ElementSize/10;
            opt.disp_on=0;
            [f,v]=ggremesh(f,v,opt);
            %check the number of connected components
            [NR, ~] = segment_connected_components(f, 'explicit');
            if NR>1
                fail = true;
            end
        end
        %convert to units of mm from units of cm
        v = v*10;
        
        if size(f,1)<Params.MaxFaces && ~fail
            valid_mesh = true;
        end
    end
end
    

if ~valid_mesh && ~Params.Keep
    f = [];
    v = [];
end

end