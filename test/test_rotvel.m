% 2018-01-22 19:30:42.923679457 +0100

dir = [1,0]

uv_earth = [1,0]
uv_sn_expect = [0,1]
uv_sn    = rotvel(dir,uv_earth)
err(1) = norm(uv_sn - uv_sn_expect)

uv_earth     = [0,1]
uv_sn_expect = [1,0]
uv_sn    = rotvel(dir,uv_earth)
err(2) = norm(uv_sn - uv_sn_expect)


dir = [0,1]

uv_earth     = [ 1,0]
uv_sn_expect = [-1,0]
uv_sn    = rotvel(dir,uv_earth)
err(3) = norm(uv_sn - uv_sn_expect)

uv_earth     = [0,1]
uv_sn_expect = [0,1]
uv_sn    = rotvel(dir,uv_earth)
err(4) = norm(uv_sn - uv_sn_expect)

if (0 ~= norm(err))
	disp('test failed')
end
