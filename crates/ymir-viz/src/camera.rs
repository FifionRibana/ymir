use bevy::input::mouse::{AccumulatedMouseMotion, AccumulatedMouseScroll};
use bevy::prelude::*;
use bevy_egui::input::EguiWantsInput;

use crate::state::CursorWorldPos;

pub struct CameraPlugin;

impl Plugin for CameraPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, setup_camera)
            .add_systems(Update, (update_cursor_world_pos, camera_pan, camera_zoom));
    }
}

#[derive(Component)]
pub struct MainCamera;

fn setup_camera(mut commands: Commands) {
    commands.spawn((Camera2d, MainCamera));
}

fn camera_pan(
    egui_input: Res<EguiWantsInput>,
    mouse: Res<ButtonInput<MouseButton>>,
    motion: Res<AccumulatedMouseMotion>,
    mut query: Query<(&mut Transform, &Projection), With<MainCamera>>,
) {
    if egui_input.is_pointer_over_area() {
        return;
    }
    if !mouse.pressed(MouseButton::Left) {
        return;
    }

    let delta = motion.delta;
    if delta == Vec2::ZERO {
        return;
    }

    if let Ok((mut transform, projection)) = query.single_mut() {
        let scale = match projection {
            Projection::Orthographic(ortho) => ortho.scale,
            _ => 1.0,
        };
        transform.translation.x -= delta.x * scale;
        transform.translation.y += delta.y * scale;
    }
}

fn camera_zoom(
    egui_input: Res<EguiWantsInput>,
    scroll: Res<AccumulatedMouseScroll>,
    mut query: Query<(&mut Transform, &mut Projection), With<MainCamera>>,
    cursor_pos: Res<CursorWorldPos>,
) {
    if egui_input.is_pointer_over_area() {
        return;
    }

    let scroll_amount = scroll.delta.y * 0.1;
    if scroll_amount == 0.0 {
        return;
    }

    let Ok((mut transform, mut projection)) = query.single_mut() else {
        return;
    };

    let Projection::Orthographic(ref mut ortho) = *projection else {
        return;
    };

    let old_scale = ortho.scale;
    let new_scale = (old_scale * (1.0 - scroll_amount)).clamp(0.05, 50.0);
    ortho.scale = new_scale;

    if let Some(cursor_world) = cursor_pos.pos {
        let factor = 1.0 - new_scale / old_scale;
        let cam_pos = transform.translation.truncate();
        let diff = cursor_world - cam_pos;
        transform.translation.x += diff.x * factor;
        transform.translation.y += diff.y * factor;
    }
}

fn update_cursor_world_pos(
    windows: Query<&Window>,
    camera_q: Query<(&Camera, &GlobalTransform), With<MainCamera>>,
    mut cursor_pos: ResMut<CursorWorldPos>,
) {
    let Ok(window) = windows.single() else {
        cursor_pos.pos = None;
        return;
    };
    let Ok((camera, camera_transform)) = camera_q.single() else {
        cursor_pos.pos = None;
        return;
    };

    cursor_pos.pos = window
        .cursor_position()
        .and_then(|cursor| camera.viewport_to_world_2d(camera_transform, cursor).ok());
}
