classdef OpenEyeSim < handle

    properties (Access = private)
        id_ % ID of the session.
    end

    methods
        function this = OpenEyeSim(filename)
            assert(ischar(filename));
            this.id_ = OpenEyeSim_('new', filename);
        end

        %DELETE Destructor.
        function delete(this)
            OpenEyeSim_('delete', this.id_);
        end

        %PUT parameters to the simulator
        function set_params(this, texture, angle, distance, strabismusAngle, planeScale)
            assert(isscalar(this));
            OpenEyeSim_('set_params', this.id_, texture, angle, distance, strabismusAngle, planeScale);
        end

        %PUT parameters to the simulator
        function add_texture(this, number, texture)
            assert(isscalar(this));
            OpenEyeSim_('add_texture', this.id_, number, texture);
        end

        %PUT Save something to the database.
        function result = generate_left(this)
            assert(isscalar(this));
            result = OpenEyeSim_('generate_left', this.id_);
        end

        %PUT Save something to the database.
        function result = generate_right(this)
            assert(isscalar(this));
            result = OpenEyeSim_('generate_right', this.id_);
        end

        %PUT Save something to the database.
        function request(this)
            assert(isscalar(this));
            OpenEyeSim_('request', this.id_);
        end

        %PUT Save something to the database.
        function initRenderer(this)
            assert(isscalar(this));
            OpenEyeSim_('initRenderer', this.id_);
        end

        %PUT Save something to the database.
        function reinitRenderer(this)
            assert(isscalar(this));
            OpenEyeSim_('reinitRenderer', this.id_);
        end
    end
end
