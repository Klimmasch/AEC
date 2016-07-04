classdef OpenEyeSimV2 < handle

properties (Access = private)
  id_ % ID of the session.
end

methods
  function this = OpenEyeSimV2(filename)  
    assert(ischar(filename));
    this.id_ = OpenEyeSimV2_('new', filename);
  end

  function delete(this)
  %DELETE Destructor.
    OpenEyeSimV2_('delete', this.id_);
  end

  function set_params(this, texture, angle, distance, strabismusAngle, planeScale)
  %PUT parameters to the simulator
    assert(isscalar(this));
    OpenEyeSimV2_('set_params', this.id_, texture, angle, distance, strabismusAngle, planeScale);
  end
  
  function add_texture(this, number, texture)
  %PUT parameters to the simulator
    assert(isscalar(this));
    OpenEyeSimV2_('add_texture', this.id_, number, texture);
  end
  
  function result = generate_left(this)
  %PUT Save something to the database.
    assert(isscalar(this));
    result = OpenEyeSimV2_('generate_left', this.id_);
  end
  
  function result = generate_right(this)
  %PUT Save something to the database.
    assert(isscalar(this));
    result = OpenEyeSimV2_('generate_right', this.id_);
  end
  
  function request(this)
  %PUT Save something to the database.
    assert(isscalar(this));
    OpenEyeSimV2_('request', this.id_);
  end
  
   function initRenderer(this)
  %PUT Save something to the database.
    assert(isscalar(this));
    OpenEyeSimV2_('initRenderer', this.id_);
   end
  
   function reinitRenderer(this)
  %PUT Save something to the database.
    assert(isscalar(this));
    OpenEyeSimV2_('reinitRenderer', this.id_);
  end
  
end

end