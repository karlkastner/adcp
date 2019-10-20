% Sun Jan  5 00:06:42 WIB 2014
% Karl Kastner, Berlin
%% process the discharge
% TODO filter in the beginning
% TODO avoid writing to the ADCP struct
% TODO error estimation is at the moment incompatible with time slicing
%      has to be performed on predicted values, not parameters,
%      matrix also becomes 3d

function [obj] = process_discharge(obj,adcp)
	% calculate the discharge
%	timer     = tic();

	switch (class(obj.grid_n))
	case {'0d'}
		obj.transect.process_discharge(adcp);
	case {'Grid1','RegularizedInterpolator1','SparseMesh1'}

	        % TODO move the T-limit here and rename into identify_t
	        obj.determine_time_slots(adcp);
	
	        % rotate the velocity vectors into cross section coordinates
	        tid   = adcp.ens.tid(:,obj.cdx);
	        msk   = logical(tid); %(tid == obj.cdx);
	        % adcp.to_cs(obj.dir,msk);
	
	        % compute depth average velocity
	        % TODO this is only required for 1d processing
		% again, this is awkward to be put here, cs should not be a member of adcp
	        % adcp.depth_average_velocity('cs');

		% generate tn-mesh, has to come before fitting of the bed-profile
	        obj.generate_mesh_tn(adcp,msk);
	
		% recomputation with partial data only necessary for error estimation
	        switch (obj.errmode)
	            case {'standard','jackknife'}
			nt = adcp.ens.n_transect;
	            otherwise
	                nt = 0;
	        end % switch
	
	        % for each transect
	        % zero must come last (all elements)
		% this loop is only iterated if error estimation is active
	        for idx=[1:nt, 0]
		
		    if (0 == idx)
			% all transects
	                ensmask = true(adcp.ens.n,1);
		    else
	                switch (lower(obj.errmode))
	                case {'standard'}
	                    % use only current transect
	                    ensmask = (adcp.ens.transect == idx);
	                case {'jackknife'}
	                    % use all but current transect
	                    % Note : the estimate for most values is the mean,
	                    %        hence the jn and standard estimate are almost the same
	                    ensmask = (adcp.ens.transect ~= idx);
	                otherwise
			    % 'none' should not come here, as loop only iterates to 0
	                    error('unimplemented error method')
	                end % switch errmode
		    end
	
	            % TODO tid should be member of adcp.ens
	            %    tmsk = (cdx == obj.tid);
	            msk = logical(adcp.ens.tid(:,obj.cdx)) & ensmask;
	            %msk = (obj.cdx == adcp.ens.tid(:,obj.cdx)) & ensmask;
	        
	            % create meshes and indices for ensembles and bins
		    % obj.build_indices(adcp,msk);
	        
	            % estimate the bottom profile of the cross section
	            % TODO regression of bottom profile for each transect individually may lead to different meshes
	            obj.fit_bed_profile(adcp, msk);
			
		    % mesh z-axis, this has to come after fitting of the bed-profile
		    obj.generate_mesh_tnz(adcp,msk);
        
	            % extrapolate the bottom profile to the unmeasured parts at the sides
		    % TODO should be process bottom profile
	            %obj.extrapolate_bed_profile();
	        
	            % use predefined bottom profile (for better treatment of side flow)
	            if (~isempty(obj.external_bottom))
	                    % keep the bottom
	                    obj.grid_n.val.bottom0 = obj.grid_n.val.bottom;
	                    % take over provided profile
	                    %obj.cs(cdx).grid_n.val.bottom = cvec(interp1(obj.cs(cdx).external_bottom.N, ...
	                    %        obj.external_bottom.val{cdx}, obj.cs(cdx).grid_n.cX1,'linear'));
	                    obj.grid_n.val.bottom = cvec(interp1(obj.external_bottom.N, ...
	                            obj.external_bottom.val, obj.grid_n.cX1,'linear'));
	            end % if
	            
	            % average velocity in mesh cells
	            obj.process_velocity_tn(adcp,msk);
	            obj.process_velocity_tnz(adcp,msk);
	        
	            % filter the velocity
	            % TODO, this should be performed on input data
	        %    for idx=1:size(obj.grid_nZ.val.U,3)
	        %        obj.grid_nZ.val.U(:,:,idx) = filter_(obj.grid_nZ.val.U(:,:,idx),5);
	        %    end
	        
	            % extrapolate bottom, top and side flow
	            %if (obj.extrapolate)
	                obj.extrapolate_velocity();
	            %end
	        
	            % calculate area and discharge of the cross section
	            % obj.integrate_discharge();
	        
	            % estimate shear velocity and roughness length
		    % TODO make own processing routine
	            obj.fit_vertical_profile_of_velocity(adcp,msk);
	        
	            % calculate auxiliary quantities
	            %obj.calc_auxiliary_quant();
	        
	            % store the intermediate results
	            if (idx>0)
	                switch (lower(obj.errmode))
	                    case {'standard','jackknife'}
	                    % vector values
	                    for field=fieldnames(obj.pseudo.grid_n)'
	                        obj.pseudo.grid_n.(field{1})(:,idx) = ...
	                                obj.grid_n.val.(field{1})(:,1);
	                    end
	                    % scalar values
	                    % TODO these values are non-scalar if tmode = function
	                    for field={'area.total', 'discharge.total', 'U', 'ln_z0'}
	                        val = getfield_deep(obj.pseudo.cs,field{1});
	                        val(end+1) = getfield_deep(obj,field{1});
	                        obj.pseudo.cs = setfield_deep(obj.pseudo.cs,field{1},val);
	                    end
	                otherwise
	                    error('unimplemented error method');
	                end % end switch errmode
	            end % if idx
	        end % for idx (transects)
	    
	        % estimate the error
	        switch (lower(obj.errmode))
	        case {'none'}
	        	% do nothing
	        case {'standard'}
	        	% vector values
	        	for field=fieldnames(obj.pseudo.grid_n)'
	                    obj.grid_n.err.(field{1}) = nanserr( obj.pseudo.grid_n.(field{1}), 2);
	                end 
	        	% scalar values
	        	for field={'area.total', 'discharge.total', 'U', 'ln_z0'}
	            		obj.err = setfield_deep(obj.err,field{1}, ...
	                	serr(getfield_deep(obj.pseudo.cs,field{1})));
	        	end
	    	case {'jackknife'}
	        	error('not jet implemented');
	    	otherwise
	        	error('unimplemented error method');
	        end % end switch errmode
	
%		runtime    = toc(timer);
%		timer      = tic();
%		runtime(2) = toc(timer);
	otherwise
		error('here');
	end % switch vmethod
end % process_discharge

