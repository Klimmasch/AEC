classdef OpenEyeSim < handle

properties (Access = private)
  id_ % ID of the session.
end

methods
  function this = OpenEyeSim(filename)  
    assert(ischar(filename));
    this.id_ = OpenEyeSim_('new', filename);
  end

  function delete(this)
  %DELETE Destructor.
    OpenEyeSim_('delete', this.id_);
  end

  function set_params(this, texture, angle, distance)
  %PUT parameters to the simulator
    assert(isscalar(this));
    OpenEyeSim_('set_params', this.id_, texture, angle, distance);
  end
  
  function add_texture(this, number, texture)
  %PUT parameters to the simulator
    assert(isscalar(this));
    OpenEyeSim_('add_texture', this.id_, number, texture);
  end
  
  function result = generate_left(this)
  %PUT Save something to the database.
    assert(isscalar(this));
    result = OpenEyeSim_('generate_left', this.id_);
  end
  
  function result = generate_right(this)
  %PUT Save something to the database.
    assert(isscalar(this));
    result = OpenEyeSim_('generate_right', this.id_);
  end
  
  function request(this)
  %PUT Save something to the database.
    assert(isscalar(this));
    OpenEyeSim_('request', this.id_);
  end
  
   function initRenderer(this)
  %PUT Save something to the database.
    assert(isscalar(this));
    OpenEyeSim_('initRenderer', this.id_);
   end
  
   function reinitRenderer(this)
  %PUT Save something to the database.
    assert(isscalar(this));
    OpenEyeSim_('reinitRenderer', this.id_);
  end
  
end

end
